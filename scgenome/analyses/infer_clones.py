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
import scgenome.breakpointdata
import scgenome.snpdata
import scgenome.snvphylo
import scgenome.db.search
import scgenome.cnclones
import scgenome.db.qc
import scgenome.loaders.qc
import scgenome.loaders.snv
import scgenome.loaders.allele
import scgenome.loaders.breakpoint

import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv
import wgs_analysis.annotation.position

import dollo
import dollo.tasks

import dbclients.tantalus
from dbclients.basicclient import NotFoundError
import datamanagement.transfer_files


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


# Threshold on total cell read counts
read_count_threshold = 500000

# SNV Calling thresholds
museq_score_threshold = None
strelka_score_threshold = None

# Cluster pruning thresholds
is_original_cluster_mean_threshold = 0.5
cluster_size_threshold = 50

# Threshold on total haplotype allele counts
# for calculating allele specific copy number
total_allele_counts_threshold = 6

# This is a global setting that should
# not be changed
cn_bin_size = 500000

results_storage_name = 'singlecellresults'


def load_cell_cycle_data(
        tantalus_api,
        jira_ticket,
        analysis_version='v0.0.2',
        results_storage_name='singlecellresults',
    ):
    """ Load cell cycle predictions for a matched hmmcopy jira ticket
    """
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    analysis = tantalus_api.get(
        'analysis',
        jira_ticket=jira_ticket,
        version=analysis_version,
        analysis_type__name='cell_state_classifier',
    )

    results = tantalus_api.get(
        'resultsdataset',
        analysis=analysis['id'],
    )

    file_instances = tantalus_api.get_dataset_file_instances(
        results['id'], 'resultsdataset', results_storage_name,
        filters={'filename__endswith': '.csv'})

    if len(file_instances) != 1:
        raise ValueError(f'found {len(file_instances)} file instances for cell cycle for {jira_ticket}')

    file_instance = file_instances[0]

    f = storage_client.open_file(file_instance['file_resource']['filename'])
    data = pd.read_csv(f)

    return data


def import_image_feature_data(
        tantalus_api,
        library_ids,
        results_storage_name='singlecellresults',
):
    """ Import image features for a list of libraries
    """
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    image_feature_data = []

    for library_id in library_ids:
        try:
            results = scgenome.db.search.search_image_feature_results(tantalus_api, library_id)
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


def retrieve_cn_data(tantalus_api, library_id, sample_ids, local_cache_directory):
    """ Retrieve comprehensive metrics data for a library
    """

    logging.info(f'library {library_id}')

    analysis = scgenome.db.search.search_hmmcopy_analysis(tantalus_api, library_id)
    jira_ticket = analysis['jira_ticket']

    scgenome.db.qc.cache_qc_results(
        analysis['jira_ticket'],
        local_cache_directory)

    ticket_directory = os.path.join(local_cache_directory, analysis['jira_ticket'])

    results = scgenome.loaders.qc.load_qc_data(ticket_directory)

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']
    align_metrics_data = results['align_metrics']

    assert isinstance(cn_data['cell_id'].dtype, pd.api.types.CategoricalDtype)
    assert isinstance(metrics_data['cell_id'].dtype, pd.api.types.CategoricalDtype)

    metrics_data = metrics_data.merge(align_metrics_data, how='left')
    assert 'cell_id_x' not in metrics_data

    if 'is_s_phase' not in metrics_data:
        cell_cycle_data = load_cell_cycle_data(
            tantalus_api,
            analysis['jira_ticket'])
        cell_cycle_data['cell_id'] = cell_cycle_data['cell_id'].astype('category')

        scgenome.utils.union_categories([
            cn_data,
            metrics_data,
            align_metrics_data,
            cell_cycle_data])

        metrics_data = metrics_data.merge(cell_cycle_data, how='left')
        assert 'cell_id_x' not in metrics_data

    assert isinstance(cn_data['cell_id'].dtype, pd.api.types.CategoricalDtype)
    assert isinstance(metrics_data['cell_id'].dtype, pd.api.types.CategoricalDtype)

    return cn_data, metrics_data


def retrieve_cn_data_multi(library_ids, sample_ids, local_cache_directory):
    tantalus_api = dbclients.tantalus.TantalusApi()

    cn_data = []
    metrics_data = []

    for library_id in library_ids:
        lib_cn_data, lib_metrics_data = retrieve_cn_data(
            tantalus_api, library_id, sample_ids, local_cache_directory)

        cn_data.append(lib_cn_data)
        metrics_data.append(lib_metrics_data)

    cn_data = scgenome.utils.concat_with_categories(cn_data)
    metrics_data = scgenome.utils.concat_with_categories(metrics_data)

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads_hmmcopy'] > read_count_threshold]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data['cell_id']
    cn_data = cn_data[cn_data['cell_id'].isin(cell_ids)]

    return cn_data, metrics_data


def retrieve_pseudobulk_data(
        ticket_id, clusters, local_cache_directory, results_prefix,
        museq_score_threshold=None, strelka_score_threshold=None,
        snvs_num_cells_threshold=2, snvs_sum_alt_threshold=2,
    ):
    """ Retrieve SNV, breakpoint and allele data
    """
    ticket_directory = os.path.join(local_cache_directory, ticket_id)

    logging.info('snv data')
    snv_results = scgenome.loaders.snv.load_snv_data(
        ticket_directory,
        museq_filter=museq_score_threshold,
        strelka_filter=strelka_score_threshold,
    )
    snv_data = snv_results['snv_data']
    snv_count_data = snv_results['snv_count_data']

    snv_data, snv_count_data = scgenome.snvdata.filter_snv_data(
        snv_data,
        snv_count_data,
        num_cells_threshold=snvs_num_cells_threshold,
        sum_alt_threshold=snvs_sum_alt_threshold,
        figures_prefix=results_prefix + 'snv_loading_',
    )

    logging.info('allele data')
    allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
        ticket_directory,
    )
    allele_data = allele_results['allele_counts']
    allele_data = scgenome.snpdata.calculate_cluster_allele_counts(
        allele_data,
        clusters,
        cn_bin_size,
    )

    logging.info('breakpoint data')
    breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data(
        ticket_directory,
    )
    breakpoint_data = breakpoint_results['breakpoint_data']
    breakpoint_count_data = breakpoint_results['breakpoint_count_data']
    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.annotate_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
    )

    logging.info('pre filter breakpoint library portrait')
    scgenome.breakpointdata.plot_library_portrait(
        breakpoint_data,
        results_prefix + 'breakpoints_unfiltered_',
    )

    logging.info('filter breakpoints')
    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.filter_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
    )

    logging.info('post filter breakpoint library portrait')
    scgenome.breakpointdata.plot_library_portrait(
        breakpoint_data,
        results_prefix + 'breakpoints_filtered_',
    )

    logging.info('plot breakpoint clustering')
    scgenome.breakpointdata.plot_breakpoint_clustering(
        breakpoint_data,
        breakpoint_count_data,
        clusters,
        results_prefix + 'clone_breakpoints_',
    )

    return snv_data, snv_count_data, allele_data, breakpoint_data, breakpoint_count_data


@click.group()
@click.pass_context
@click.argument('results_prefix')
@click.argument('local_cache_directory')
def infer_clones_cmd(ctx, results_prefix, local_cache_directory):
    ctx.obj['results_prefix'] = results_prefix
    ctx.obj['local_cache_directory'] = local_cache_directory


@infer_clones_cmd.command('retriever-cn')
@click.pass_context
@click.option('--library_id')
@click.option('--sample_id')
@click.option('--library_ids_filename')
@click.option('--sample_ids_filename')
def retrieve_cn_cmd(
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

    retrieve_cn(library_ids, sample_ids, results_prefix, local_cache_directory)


def retrieve_cn(library_ids, sample_ids, results_prefix, local_cache_directory):
    logging.info('retrieving cn data')
    cn_data, metrics_data = retrieve_cn_data_multi(
        library_ids,
        sample_ids,
        local_cache_directory,
    )

    cn_data.to_pickle(results_prefix + 'cn_data.pickle')
    metrics_data.to_pickle(results_prefix + 'metrics_data.pickle')


@infer_clones_cmd.command('cluster-cn')
@click.pass_context
def cluster_cn_cmd(ctx):
    results_prefix = ctx.obj['results_prefix']
    cluster_cn(results_prefix)


def cluster_cn(results_prefix, cluster_size_threshold=50):
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
        cluster_size_threshold=cluster_size_threshold,
    )

    mitotic_errors = scgenome.cnclones.calculate_mitotic_errors(
        cn_data,
        clusters,
        results_prefix + 'mitotic_errors_',
    )

    clusters.to_pickle(results_prefix + 'clusters.pickle')
    filter_metrics.to_pickle(results_prefix + 'filter_metrics.pickle')
    cell_clone_distances.to_pickle(results_prefix + 'cell_clone_distances.pickle')
    final_clusters.to_pickle(results_prefix + 'final_clusters.pickle')
    mitotic_errors.to_pickle(results_prefix + 'mitotic_errors.pickle')


@infer_clones_cmd.command('pseudobulk-analysis')
@click.pass_context
@click.argument('pseudobulk_ticket')
def pseudobulk_analysis_cmd(ctx, pseudobulk_ticket):
    results_prefix = ctx.obj['results_prefix']
    local_cache_directory = ctx.obj['local_cache_directory']
    pseudobulk_analysis(pseudobulk_ticket, results_prefix, local_cache_directory)


def pseudobulk_analysis(pseudobulk_ticket, results_prefix, local_cache_directory):
    cn_data = pd.read_pickle(results_prefix + 'cn_data.pickle')
    clusters = pd.read_pickle(results_prefix + 'clusters.pickle')
    final_clusters = pd.read_pickle(results_prefix + 'final_clusters.pickle')

    logging.info('retrieving pseudobulk data')
    snv_data, snv_count_data, allele_data, breakpoint_data, breakpoint_count_data = retrieve_pseudobulk_data(
        pseudobulk_ticket,
        final_clusters,
        local_cache_directory,
        results_prefix + 'retrieve_pseudobulk_data_',
    )

    logging.info('calculate cluster allele cn')
    allele_cn = scgenome.snpdata.calculate_cluster_allele_cn(
        cn_data,
        allele_data,
        clusters,
        results_prefix + 'calculate_cluster_allele_cn_',
    )
    
    logging.info('bulk snv analysis')
    scgenome.snvdata.run_bulk_snv_analysis(
        snv_data,
        snv_count_data,
        final_clusters[['cell_id']].drop_duplicates(),
        results_prefix + 'run_bulk_snv_analysis_',
    )

    logging.info('snv phylogenetics')
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
