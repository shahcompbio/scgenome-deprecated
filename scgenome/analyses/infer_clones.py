import sys
import os
import logging
import click
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functools
import itertools
import pickle

import seaborn
import numpy as np
import pandas as pd
import pylab
import sklearn.preprocessing
import scipy.spatial.distance

import scgenome
import scgenome.utils
import scgenome.cncluster
import scgenome.cnplot
import scgenome.snvdata
import scgenome.snpdata
import scgenome.snvphylo
import scgenome.dbsearch
import scgenome.hmmcopy
import scgenome.cnclones
import scgenome.pseudobulk

import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv
import wgs_analysis.annotation.position

import dollo
import dollo.tasks

import dbclients.tantalus
from dbclients.basicclient import NotFoundError
import datamanagement.transfer_files


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


# SNV Calling thresholds
museq_score_threshold = None
strelka_score_threshold = None
snvs_num_cells_threshold = 2
snvs_sum_alt_threshold = 2

# Cluster pruning thresholds
is_original_cluster_mean_threshold = 0.5
cluster_size_threshold = 50

# Threshold on total haplotype allele counts
# for calculating allele specific copy number
total_allele_counts_threshold = 6

# This is a global setting that should
# not be changed
cn_bin_size = 500000

results_storage_name = 'singlecellblob_results'


def search_hmmcopy_analyses(
        tantalus_api,
        library_ids,
        aligner_name='BWA_ALN_0_5_7',
):
    """ Search for hmmcopy results and analyses for a list of libraries.
    """
    hmmcopy_results = {}
    hmmcopy_tickets = {}

    for library_id in library_ids:
        results, analysis = scgenome.dbsearch.search_hmmcopy_results(
            tantalus_api, library_id, aligner_name=aligner_name)
        hmmcopy_results[library_id] = results
        hmmcopy_tickets[library_id] = analysis['jira_ticket']

    return hmmcopy_results, hmmcopy_tickets


def import_cell_cycle_data(
        tantalus_api,
        library_ids,
        hmmcopy_results,
        results_storage_name='singlecellblob_results',
):
    """ Import cell cycle predictions for a list of libraries
    """
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    cell_cycle_data = []

    for library_id in library_ids:
        results, analysis = scgenome.dbsearch.search_cell_cycle_results(
            tantalus_api, library_id, hmmcopy_results[library_id])

        file_instances = tantalus_api.get_dataset_file_instances(
            results['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f)
            data['library_id'] = library_id
            cell_cycle_data.append(data)

    cell_cycle_data = pd.concat(cell_cycle_data, ignore_index=True, sort=True)

    return cell_cycle_data


def import_image_feature_data(
        tantalus_api,
        library_ids,
        results_storage_name='singlecellblob_results',
):
    """ Import image features for a list of libraries
    """
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    image_feature_data = []

    for library_id in library_ids:
        try:
            results = scgenome.dbsearch.search_image_feature_results(tantalus_api, library_id)
        except NotFoundError:
            logging.info('no image data for {}'.format(library_id))
            continue
        file_instances = tantalus_api.get_dataset_file_instances(
            results['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f, index_col=0)
            data['library_id'] = library_id
            image_feature_data.append(data)

    if len(image_feature_data) == 0:
        return pd.DataFrame()

    image_feature_data = pd.concat(image_feature_data, ignore_index=True, sort=True)

    return image_feature_data


def retrieve_cn_data(library_ids, sample_ids, local_cache_directory, results_prefix):
    tantalus_api = dbclients.tantalus.TantalusApi()

    hmmcopy_results, hmmcopy_tickets = search_hmmcopy_analyses(tantalus_api, library_ids)

    cn_data = []
    metrics_data = []

    for library_id in library_ids:
        hmmcopy = scgenome.hmmcopy.HMMCopyData(hmmcopy_tickets[library_id], local_cache_directory)
        hmmcopy_data = hmmcopy.load_cn_data(sample_ids=sample_ids)

        cn_data.append(hmmcopy_data['hmmcopy_reads'])
        metrics_data.append(hmmcopy_data['hmmcopy_metrics'])

    cn_data = scgenome.utils.concat_with_categories(cn_data)
    metrics_data = scgenome.utils.concat_with_categories(metrics_data)

    cell_cycle_data = import_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)
    cell_cycle_data['cell_id'] = pd.Categorical(cell_cycle_data['cell_id'], categories=metrics_data['cell_id'].cat.categories)
    metrics_data = metrics_data.merge(cell_cycle_data)

    image_feature_data = import_image_feature_data(tantalus_api, library_ids)

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads_hmmcopy'] > 500000]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data['cell_id']
    cn_data = cn_data[cn_data['cell_id'].isin(cell_ids)]

    # TODO: Remove temporary fixup
    if 'total_mapped_reads_hmmcopy' not in metrics_data:
         metrics_data['total_mapped_reads_hmmcopy'] = metrics_data['total_mapped_reads']
    elif metrics_data['total_mapped_reads_hmmcopy'].isnull().any():
        fix_read_count = metrics_data['total_mapped_reads_hmmcopy'].isnull()
        metrics_data.loc[fix_read_count, 'total_mapped_reads_hmmcopy'] = (
            metrics_data.loc[fix_read_count, 'total_mapped_reads'])

    return cn_data, metrics_data, image_feature_data


def retrieve_pseudobulk_data(
        ticket_id, clusters, local_cache_directory, results_prefix,
        museq_score_threshold=None, strelka_score_threshold=None,
    ):
    """ Retrieve SNV, breakpoint and allele data
    """

    pseudobulk = scgenome.pseudobulk.PseudobulkData(ticket_id, local_cache_directory)

    snv_data, snv_count_data = scgenome.snvdata.load_snv_data(
        pseudobulk,
        museq_filter=museq_score_threshold,
        strelka_filter=strelka_score_threshold,
        num_cells_threshold=snvs_num_cells_threshold,
        sum_alt_threshold=snvs_sum_alt_threshold,
        figures_prefix=results_prefix + 'snv_loading_',
    )

    allele_data = pseudobulk.load_haplotype_allele_counts()
    allele_data = scgenome.snpdata.calculate_cluster_allele_counts(allele_data, clusters, cn_bin_size)

    breakpoint_data, breakpoint_count_data = pseudobulk.load_breakpoint_data()

    return snv_data, snv_count_data, allele_data, breakpoint_data, breakpoint_count_data


@click.group()
@click.pass_context
@click.argument('results_prefix')
@click.argument('local_cache_directory')
def infer_clones_cmd(ctx, results_prefix, local_cache_directory):
    ctx.obj['results_prefix'] = results_prefix
    ctx.obj['local_cache_directory'] = local_cache_directory


@infer_clones_cmd.command()
@click.pass_context
@click.option('--library_id')
@click.option('--sample_id')
@click.option('--library_ids_filename')
@click.option('--sample_ids_filename')
def retrieve_cn(
        ctx,
        library_id=None,
        library_ids_filename=None,
        sample_id=None,
        sample_ids_filename=None,
    ):

    results_prefix = ctx.obj['results_prefix']
    local_cache_directory = ctx.obj['local_cache_directory']

    if library_id is not None:
        library_ids = [library_id]
    elif library_ids_filename is not None:
        library_ids = [l.strip() for l in open(library_ids_filename).readlines()]
    else:
        raise Exception('must specify library_id or library_ids_filename')

    if sample_id is not None:
        sample_ids = [sample_id]
    elif sample_ids_filename is not None:
        sample_ids = [l.strip() for l in open(sample_ids_filename).readlines()]
    else:
        raise Exception('must specify sample_id or sample_ids_filename')

    logging.info('retrieving cn data')
    cn_data, metrics_data, image_feature_data = retrieve_cn_data(
        library_ids,
        sample_ids,
        local_cache_directory,
        results_prefix + 'retrieve_cn_',
    )

    cn_data.to_pickle(results_prefix + 'cn_data.pickle')
    metrics_data.to_pickle(results_prefix + 'metrics_data.pickle')
    image_feature_data.to_pickle(results_prefix + 'image_feature_data.pickle')


@infer_clones_cmd.command()
@click.pass_context
def cluster_cn(ctx):

    results_prefix = ctx.obj['results_prefix']

    cn_data = pd.read_pickle(results_prefix + 'cn_data.pickle')
    metrics_data = pd.read_pickle(results_prefix + 'metrics_data.pickle')

    logging.info('calculating clusters')
    clusters, filter_metrics = scgenome.cnclones.calculate_clusters(
        cn_data,
        metrics_data,
        results_prefix + 'calculate_clusters_',
    )

    cell_clone_distances = scgenome.cnclones.calculate_cell_clone_distances(
        cn_data,
        clusters,
        results_prefix + 'calculate_cell_clone_distances_',
    )

    final_clusters = scgenome.cnclones.finalize_clusters(
        cn_data,
        metrics_data,
        clusters,
        filter_metrics,
        cell_clone_distances,
        results_prefix + 'finalize_clusters_',
    )

    clusters.to_pickle(results_prefix + 'clusters.pickle')
    filter_metrics.to_pickle(results_prefix + 'filter_metrics.pickle')
    cell_clone_distances.to_pickle(results_prefix + 'cell_clone_distances.pickle')
    final_clusters.to_pickle(results_prefix + 'final_clusters.pickle')


@infer_clones_cmd.command()
@click.pass_context
@click.argument('pseudobulk_ticket')
def pseudobulk_analysis(ctx, pseudobulk_ticket):

    results_prefix = ctx.obj['results_prefix']
    local_cache_directory = ctx.obj['local_cache_directory']

    cn_data = pd.read_pickle(results_prefix + 'cn_data.pickle')
    clusters = pd.read_pickle(results_prefix + 'clusters.pickle')
    final_clusters = pd.read_pickle(results_prefix + 'final_clusters.pickle')

    snv_data, snv_count_data, allele_data, breakpoint_data, breakpoint_count_data = retrieve_pseudobulk_data(
        pseudobulk_ticket,
        final_clusters,
        local_cache_directory,
        results_prefix + 'retrieve_pseudobulk_data_',
    )

    allele_cn = scgenome.snpdata.calculate_cluster_allele_cn(
        cn_data,
        allele_data,
        clusters,
        results_prefix + 'calculate_cluster_allele_cn_',
    )
    
    scgenome.snvdata.run_bulk_snv_analysis(
        snv_data,
        snv_count_data,
        final_clusters[['cell_id']].drop_duplicates(),
        results_prefix + 'run_bulk_snv_analysis_',
    )

    snv_ml_tree, snv_tree_annotations = scgenome.snvphylo.run_snv_phylogenetics(
        snv_count_data,
        allele_cn,
        final_clusters,
        results_prefix + 'run_snv_phylogenetics_', 
    )

    snv_data.to_pickle(results_prefix + 'snv_data.pickle')
    snv_count_data.to_pickle(results_prefix + 'snv_count_data.pickle')
    allele_data.to_pickle(results_prefix + 'allele_data.pickle')
    breakpoint_data.to_pickle(results_prefix + 'breakpoint_data.pickle')
    breakpoint_count_data.to_pickle(results_prefix + 'breakpoint_count_data.pickle')

    allele_cn.to_pickle(results_prefix + 'breakpoint_data.pickle')
    with open(results_prefix + 'snv_ml_tree.pickle', 'wb') as f:
        pickle.dump(snv_ml_tree, f)
    snv_tree_annotations.to_pickle(results_prefix + 'snv_tree_annotations.pickle')


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    infer_clones_cmd(obj={})
