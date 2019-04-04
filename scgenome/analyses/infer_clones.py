import sys
import os
import logging
import click
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wget
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
import scgenome.dataimport
import scgenome.cncluster
import scgenome.cnplot

import dbclients.tantalus
from dbclients.basicclient import NotFoundError
import datamanagement.transfer_files


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


results_storage_name = 'singlecellblob_results'


def get_hmmcopy_analyses(tantalus_api, library_ids):
    hmmcopy_results = {}
    hmmcopy_tickets = []

    for library_id in library_ids:
        logging.info('hmmcopy data for {}'.format(library_id))

        analyses = list(tantalus_api.list(
            'analysis',
            analysis_type__name='hmmcopy',
            input_datasets__library__library_id=library_id,
        ))
        aln_analyses = []
        for analysis in analyses:
            aligners = set()
            is_complete = True
            for dataset_id in analysis['input_datasets']:
                dataset = tantalus_api.get('sequencedataset', id=dataset_id)
                if not dataset['is_complete']:
                    is_complete = False
                aligners.add(dataset['aligner'])
            if len(aligners) != 1:
                # HACK: should be an exception but will remove when datasets are cleaned up
                continue
                raise Exception('found {} aligners for analysis {}'.format(
                    len(aligners), analysis['id']))
            aligner = aligners.pop()
            if is_complete and aligner == 'BWA_ALN_0_5_7':
                aln_analyses.append(analysis)
        if len(aln_analyses) != 1:
            raise Exception('found {} hmmcopy analyses for {}: {}'.format(
                len(aln_analyses), library_id, [a['id'] for a in aln_analyses]))
        analysis = aln_analyses[0]
        results = tantalus_api.get(
            'resultsdataset',
            analysis=analysis['id'],
        )
        hmmcopy_results[library_id] = results
        hmmcopy_tickets.append(analysis['jira_ticket'])

    return hmmcopy_results, hmmcopy_tickets


def get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results):
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    cell_cycle_data = []

    for library_id in library_ids:
        logging.info('cell cycle data for {}'.format(library_id))

        classifier = tantalus_api.get(
            'analysis',
            analysis_type='cell_state_classifier',
            version='v0.0.1',
            input_results__id=hmmcopy_results[library_id]['id'],
        )
        features = tantalus_api.get(
            'resultsdataset',
            analysis=classifier['id'],
        )
        file_instances = tantalus_api.get_dataset_file_instances(
            features['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f)
            data['library_id'] = library_id
            cell_cycle_data.append(data)
    cell_cycle_data = pd.concat(cell_cycle_data, ignore_index=True, sort=True)

    return cell_cycle_data


def get_image_feature_data(tantalus_api, library_ids):
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    image_feature_data = []

    for library_id in library_ids:
        try:
            features = tantalus_api.get(
                'results',
                results_type='CELLENONE_FEATURES',
                results_version='v0.0.1',
                libraries__library_id=library_id,
            )
        except NotFoundError:
            logging.info('no image data for {}'.format(library_id))
            continue
        file_instances = tantalus_api.get_dataset_file_instances(
            features['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f, index_col=0)
            data['library_id'] = library_id
            image_feature_data.append(data)
    image_feature_data = pd.concat(image_feature_data, ignore_index=True, sort=True)

    return image_feature_data


def retrieve_data(tantalus_api, library_ids, sample_ids, local_storage_directory, results_prefix):
    hmmcopy_results, hmmcopy_tickets = get_hmmcopy_analyses(tantalus_api, library_ids)

    results = scgenome.dataimport.import_cn_data(
        hmmcopy_tickets,
        local_storage_directory,
        sample_ids=sample_ids,
        #ploidy_solution='2',
        #subsample=0.25,
    )

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']

    cell_cycle_data = get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)
    cell_cycle_data['cell_id'] = pd.Categorical(cell_cycle_data['cell_id'], categories=metrics_data['cell_id'].cat.categories)
    metrics_data = metrics_data.merge(cell_cycle_data)

    image_feature_data = get_image_feature_data(tantalus_api, library_ids)

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads_hmmcopy'] > 500000]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data['cell_id']
    cn_data = cn_data[cn_data['cell_id'].isin(cell_ids)]

    return cn_data, metrics_data, image_feature_data


def infer_clones(cn_data, metrics_data, results_prefix):
    """ Infer clones from copy number data.
    """

    metrics_data['filter_quality'] = (metrics_data['quality'] > 0.75)
    metrics_data['filter_reads'] = (metrics_data['total_mapped_reads_hmmcopy'] > 500000)

    # Calculate separation between predicted and normalized copy number
    cn_data['copy_state_diff'] = np.absolute(cn_data['copy'] - cn_data['state'])
    copy_state_diff = (cn_data[['cell_id', 'copy_state_diff']]
        .dropna().groupby('cell_id')['copy_state_diff']
        .mean().reset_index().dropna())
    metrics_data = metrics_data.merge(copy_state_diff)
    metrics_data['filter_copy_state_diff'] = (metrics_data['copy_state_diff'] < 1.)

    # Remove s phase cells
    # Remove low quality cells
    # Remove low coverage cells
    # Remove cells with a large divergence between copy state and norm copy number
    filtered_cells = metrics_data.loc[
        (~metrics_data['is_s_phase']) &
        metrics_data['filter_quality'] &
        metrics_data['filter_reads'] &
        metrics_data['filter_copy_state_diff'],
        ['cell_id']]

    logging.info('filtering {} of {} cells'.format(
        len(filtered_cells.index), len(metrics.index)))
    cn_data = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())
    assert isinstance(cn_data['cell_id'].dtype, pd.api.types.CategoricalDtype)

    logging.info('creating copy number matrix')
    cn = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0)
    )

    logging.info('clustering copy number')
    clusters = scgenome.cncluster.umap_hdbscan_cluster(cn)

    logging.info('merging clusters')
    cn_data = cn_data.merge(clusters[['cell_id', 'cluster_id']].drop_duplicates())

    logging.info('plotting clusters to {}*'.format(results_prefix + '_initial'))
    plot_clones(cn_data, 'cluster_id', results_prefix + '_initial')

    filter_metrics = metrics_data[[
        'cell_id',
        'is_s_phase',
        'copy_state_diff',
        'filter_quality',
        'filter_reads',
        'filter_copy_state_diff',
    ]]

    return clusters, filter_metrics


def filter_clusters(metrics_data, clusters, filter_metrics, cell_clone_distances, results_prefix):
    """ Filter clusters
    """

    filter_metrics = filter_metrics.merge(cell_clone_distances, on='cell_id')

    # Calculate whether the clusters are stable when using correlation
    # vs cityblock distance as a distance metric.
    # Note: this will filter cells that have been either erroneously or
    # correctly assigned the correct ploidy, and have thus been erroneously
    # assigned to the wrong cluster.  It will successfully filter cells that
    # have a different ploidy than there actual clone.
    filter_metrics['filter_same_cluster'] = (
        (filter_metrics['cluster_id_pearsonr'] == filter_metrics['cluster_id_cityblock'])
    )

    # Strict initial filtering, excluding s phase cells
    filter_metrics['filter_final'] = (
        ~(filter_metrics['is_s_phase']) &
        filter_metrics['filter_quality'] &
        filter_metrics['filter_reads'] &
        filter_metrics['filter_same_cluster'] &
        filter_metrics['filter_copy_state_diff'])

    # Add back in 
    filter_metrics['filter_final'] |= (
        filter_metrics['is_s_phase'] &
        filter_metrics['filter_reads'] &
        filter_metrics['filter_copy_state_diff'])

    return filter_metrics


def recalculate_distances(distance_metric, distance_method, clone_cn_matrix, cell_cn_matrix):
    """ Recalculate distances to closest cluster using some metric.
    """

    logging.info('Calculating clone cell {} distance'.format(distance_metric))
    cell_clone_corr = {}
    for cluster_id in clone_cn_matrix.columns:
        logging.info('Calculating distance for clone {}'.format(cluster_id))
        cell_clone_corr[cluster_id] = cell_cn_matrix.corrwith(
            clone_cn_matrix[cluster_id], method=distance_method)

    distance = pd.DataFrame(cell_clone_corr)
    distance.columns.name = 'cluster_id_' + distance_metric

    reclusters = pd.DataFrame(cell_clone_corr).idxmin(axis=1).dropna().astype(int)
    reclusters.name = 'cluster_id_' + distance_metric
    reclusters = reclusters.reset_index()

    reclusters = reclusters.merge(distance.stack().rename('distance_' + distance_metric).reset_index())

    return reclusters


def calculate_cell_clone_distances(cn_data, clusters, results_prefix):
    """ Calculate the distance to the closest clone for multiple metrics.
    """

    logging.info('Create clone copy number table')
    clone_cn_data = (
        cn_data
            .merge(clusters[['cell_id', 'cluster_id']].drop_duplicates())
            .groupby(['chr', 'start', 'end', 'cluster_id'])
            .agg({'copy': np.mean, 'state': np.median})
            .reset_index()
    )
    clone_cn_data['state'] = clone_cn_data['state'].round().astype(int)

    logging.info('Create matrix of cn data for all cells')
    cell_cn_matrix = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level=['cell_id']).fillna(0)
    )

    logging.info('Create a matrix of cn data for filtered clones')
    clone_cn_matrix = (
        clone_cn_data
            .set_index(['chr', 'start', 'end', 'cluster_id'])['copy']
            .unstack(level=['cluster_id']).fillna(0)
    )

    pearsonr_distances = recalculate_distances(
        'pearsonr',
        lambda u, v: 1. - scipy.stats.pearsonr(u, v)[0],
        clone_cn_matrix,
        cell_cn_matrix,
    )

    spearmanr_distances = recalculate_distances(
        'spearmanr',
        lambda u, v: 1. - scipy.stats.spearmanr(u, v)[0],
        clone_cn_matrix,
        cell_cn_matrix,
    )

    cityblock_distances = recalculate_distances(
        'cityblock',
        scipy.spatial.distance.cityblock,
        clone_cn_matrix,
        cell_cn_matrix,
    )

    clone_cell_distances = pd.concat([
        pearsonr_distances.set_index('cell_id'),
        spearmanr_distances.set_index('cell_id'),
        cityblock_distances.set_index('cell_id'),
    ], axis=1).reset_index()

    return clone_cell_distances


def plot_clones(cn_data, cluster_col, plots_prefix):
    plot_data = cn_data.copy()
    bin_filter = (plot_data['gc'] <= 0) | (plot_data['copy'].isnull())
    plot_data.loc[bin_filter, 'state'] = 0
    plot_data.loc[plot_data['copy'] > 5, 'copy'] = 5.
    plot_data.loc[plot_data['copy'] < 0, 'copy'] = 0.

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'copy', cluster_field_name=cluster_col, raw=True)
    fig.savefig(plots_prefix + '_raw_cn.pdf')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    fig.savefig(plots_prefix + '_cn_state.pdf')


def memoize(cache_filename, func, *args, **kwargs):
    if os.path.exists(cache_filename):
        logging.info('reading existing data from {}'.format(cache_filename))
        with open(cache_filename, 'rb') as f:
            data = pickle.load(f)
    else:
        data = func(*args, **kwargs)
        logging.info('writing data to {}'.format(cache_filename))
        with open(cache_filename, 'wb') as f:
            pickle.dump(data, f, protocol=4)
    return data


@click.command()
@click.argument('library_ids_filename')
@click.argument('sample_ids_filename')
@click.argument('results_prefix')
@click.argument('local_storage_directory')
def infer_clones_cmd(library_ids_filename, sample_ids_filename, results_prefix, local_storage_directory):
    tantalus_api = dbclients.tantalus.TantalusApi()

    library_ids = [l.strip() for l in open(library_ids_filename).readlines()]
    sample_ids = [l.strip() for l in open(sample_ids_filename).readlines()]

    raw_data_filename = results_prefix + '_raw_data.pickle'

    logging.info('retrieving cn data')
    cn_data, metrics_data, image_feature_data = memoize(
        raw_data_filename,
        retrieve_data,
        tantalus_api,
        library_ids,
        sample_ids,
        local_storage_directory,
        results_prefix,
    )

    clones_filename = results_prefix + '_clones.pickle'
    logging.info('inferring clones')
    shape_check = cn_data.shape
    logging.info('cn_data shape {}'.format(shape_check))
    clusters, filter_metrics = memoize(
        clones_filename,
        infer_clones,
        cn_data,
        metrics_data,
        results_prefix,
    )
    assert cn_data.shape == shape_check

    cell_clone_distances_filename = results_prefix + '_cell_clone_distances.pickle'
    cell_clone_distances = memoize(
        cell_clone_distances_filename,
        calculate_cell_clone_distances,
        cn_data,
        clusters,
        results_prefix,
    )

    filtered_clusters_filename = results_prefix + '_filtered_clones.pickle'
    filtered_clusters = memoize(
        filtered_clusters_filename,
        filter_clusters,
        metrics_data,
        clusters,
        filter_metrics,
        cell_clone_distances,
        results_prefix,
    )

    print(cn_data.head())
    print(metrics_data.head())
    print(image_feature_data.head())
    print(clone_cn_data.head())
    print(clusters.head())


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones_cmd()
