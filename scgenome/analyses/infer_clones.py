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


def retrieve_cn_data(tantalus_api, library_ids, sample_ids, local_storage_directory, results_prefix):
    hmmcopy_results, hmmcopy_tickets = scgenome.dataimport.search_hmmcopy_analyses(tantalus_api, library_ids)

    results = scgenome.dataimport.import_cn_data(
        hmmcopy_tickets,
        local_storage_directory,
        sample_ids=sample_ids,
    )

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']

    cell_cycle_data = scgenome.dataimport.import_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)
    cell_cycle_data['cell_id'] = pd.Categorical(cell_cycle_data['cell_id'], categories=metrics_data['cell_id'].cat.categories)
    metrics_data = metrics_data.merge(cell_cycle_data)

    image_feature_data = scgenome.dataimport.import_image_feature_data(tantalus_api, library_ids)

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads_hmmcopy'] > 500000]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data['cell_id']
    cn_data = cn_data[cn_data['cell_id'].isin(cell_ids)]

    return cn_data, metrics_data, image_feature_data


def calc_prop_hom_del(states):
    cndist = states.value_counts()
    cndist = cndist / cndist.sum()
    if 0 not in cndist:
        return 0
    return cndist[0]


def calculate_clusters(cn_data, metrics_data, results_prefix):
    """ Cluster copy number data.
    """

    metrics_data['filter_quality'] = (metrics_data['quality'] > 0.75)
    metrics_data['filter_reads'] = (metrics_data['total_mapped_reads_hmmcopy'] > 500000)

    # Calculate proportion homozygous deletion state
    prop_hom_del = cn_data.groupby('cell_id')['state'].apply(calc_prop_hom_del).rename('prop_hom_del').reset_index()
    metrics_data = metrics_data.merge(prop_hom_del, how='left')
    metrics_data['prop_hom_del'] = metrics_data['prop_hom_del'].fillna(0)
    metrics_data['zscore_prop_hom_del'] = scipy.stats.zscore(metrics_data['prop_hom_del'])
    metrics_data['filter_prop_hom_del'] = (metrics_data['zscore_prop_hom_del'] < 3.)

    # Calculate separation between predicted and normalized copy number
    copy_state_diff = cn_data[['cell_id', 'copy', 'state']].copy()
    copy_state_diff['copy_state_diff'] = np.absolute(copy_state_diff['copy'] - copy_state_diff['state'])
    copy_state_diff = (copy_state_diff[['cell_id', 'copy_state_diff']]
        .dropna().groupby('cell_id')['copy_state_diff']
        .mean().reset_index().dropna())
    metrics_data = metrics_data.merge(copy_state_diff)
    metrics_data['filter_copy_state_diff'] = (metrics_data['copy_state_diff'] < 1.)

    # Remove s phase cells
    # Remove low quality cells
    # Remove low coverage cells
    # Remove cells with a large divergence between copy state and norm copy number
    # Remove cells with outlier proportion of homozygous deletion
    filtered_cells = metrics_data.loc[
        (~metrics_data['is_s_phase']) &
        metrics_data['filter_quality'] &
        metrics_data['filter_reads'] &
        metrics_data['filter_copy_state_diff'] &
        metrics_data['filter_prop_hom_del'],
        ['cell_id']]

    logging.info('filtering {} of {} cells'.format(
        len(filtered_cells.index), len(metrics_data.index)))
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
        'prop_hom_del',
        'zscore_prop_hom_del',
        'filter_prop_hom_del',
    ]]

    return clusters, filter_metrics


def breakpoint_filter(metrics_data, clusters, max_breakpoints, results_prefix):
    """ Filter clusters based on breakpoint counts
    """
    breakpoint_data = (
        metrics_data
        .merge(clusters[['cell_id', 'cluster_id']])
        .groupby('cluster_id')['breakpoints']
        .mean().reset_index()
        .sort_values('breakpoints'))

    fig = plt.figure(figsize=(4, 4))
    breakpoint_data['breakpoints'].hist(bins=40)
    fig.savefig(results_prefix + '_breakpoint_hist.pdf', bbox_inches='tight')

    breakpoint_data = breakpoint_data[breakpoint_data['breakpoints'] <= max_breakpoints]

    clusters = clusters[clusters['cluster_id'].isin(breakpoint_data['cluster_id'])]

    return clusters


is_original_cluster_mean_threshold = 0.5
cluster_size_threshold = 50

def finalize_clusters(cn_data, metrics_data, clusters, filter_metrics, cell_clone_distances, results_prefix):
    """ Generate finalized filtered clusters
    """

    # Calculate the cluster assignment based on correlation
    correlation_metric = 'pearsonr'
    correlation_cluster = (
        cell_clone_distances
        .set_index(['cluster_id', 'cell_id'])[correlation_metric]
        .unstack().idxmin().rename(correlation_metric + '_cluster_id').reset_index())

    # Calculate which cells cluster assignment matches highest correlation cluster
    cluster_annotation = cell_clone_distances.merge(clusters[['cell_id', 'cluster_id']])
    cluster_annotation = cluster_annotation.merge(correlation_cluster)
    cluster_annotation['is_original'] = (cluster_annotation['cluster_id'] == cluster_annotation[correlation_metric + '_cluster_id'])

    # Plot the cityblock distance distribution of each cluster separated
    # by whether the assigned cluster equals the correlation cluster
    plot_metric = 'cityblock'
    g = seaborn.factorplot(
        x='cluster_id', y=plot_metric,
        hue='is_original', kind='strip',
        dodge=True, data=cluster_annotation, aspect=3)
    g.fig.savefig(results_prefix + '_cluster_cityblock_distance.pdf', bbox_inches='tight')

    # Calculate the proportion of each cluster that would be assigned to
    # that cluster by maximizing correlation
    cluster_annotation['is_original_f'] = cluster_annotation['is_original'] * 1.
    is_original_mean = cluster_annotation.groupby('cluster_id')['is_original_f'].mean().rename('is_original_cluster_mean').reset_index()
    cluster_annotation = cluster_annotation.merge(is_original_mean)

    # Filter cells that are not assigned to the same cluster they
    # are most correlated with
    cluster_annotation = cluster_annotation.query('is_original')

    # Filter clusters for which more than a given proportion of cells are
    # assigned to a different cluster than that which they are most
    # correlated to
    cluster_annotation = cluster_annotation.query(
        'is_original_cluster_mean > {}'.format(is_original_cluster_mean_threshold))

    # Filter clusters smaller than a given size
    cluster_annotation = cluster_annotation.merge(
        cluster_annotation.groupby('cluster_id').size().rename('cluster_size').reset_index())
    cluster_annotation = cluster_annotation.query(
        'cluster_size >= {}'.format(cluster_size_threshold))

    # Assign s phase cells to the cluster they are most correlated with
    cell_filtered_clone_distances = cell_clone_distances.merge(
        cluster_annotation[['cluster_id']].drop_duplicates())
    s_phase_cluster = cell_filtered_clone_distances.merge(
        metrics_data.query('is_s_phase')[['cell_id']])
    s_phase_cluster = (
        s_phase_cluster
        .set_index(['cell_id', 'cluster_id'])['pearsonr']
        .unstack(level=['cluster_id']).idxmin(axis=1).rename('cluster_id').reset_index())

    # Filter s phase cells
    s_phase_filter = (filter_metrics
        .query('is_s_phase')
        .query('filter_reads')
        .query('filter_copy_state_diff')
        .query('filter_prop_hom_del'))[['cell_id']]
    s_phase_cluster = s_phase_cluster.merge(s_phase_filter)

    # Create a merged set of cluster calls
    final_clusters = pd.concat([
        cluster_annotation[['cell_id', 'cluster_id']],
        s_phase_cluster[['cell_id', 'cluster_id']],
    ], ignore_index=True)
    assert not final_clusters['cell_id'].duplicated().any()

    # Plotting
    #
    
    # Plot final clusters heatmap
    logging.info('plotting clusters to {}*'.format(results_prefix + '_filter_final'))
    plot_cn_data = cn_data.merge(
        final_clusters[['cell_id', 'cluster_id']])
    plot_clones(plot_cn_data, 'cluster_id', results_prefix + '_filter_final')

    # Plot s phase proportions
    s_plot_data = (
        metrics_data
        .merge(final_clusters[['cell_id', 'cluster_id']].drop_duplicates())
        .groupby('cluster_id').agg({'is_s_phase': (np.sum, len, np.mean)}).reset_index())
    s_plot_data.columns = ['clone', 'sum', 'len', 'proportion']

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(211)
    seaborn.barplot(ax=ax, x='clone', y='proportion', data=s_plot_data, color='0.5')
    seaborn.despine()
    ax = fig.add_subplot(212)
    seaborn.barplot(ax=ax, x='clone', y='len', data=s_plot_data, color='0.5')
    seaborn.despine()
    plt.tight_layout()
    fig.savefig(results_prefix + '_clone_s_phase.pdf', bbox_inches='tight')

    return final_clusters


def recalculate_distances(distance_metric, distance_method, clone_cn_matrix, cell_cn_matrix):
    """ Recalculate distances to closest cluster using some metric.
    """

    logging.info('Calculating clone cell {} distance'.format(distance_metric))
    cell_clone_corr = {}
    for cluster_id in clone_cn_matrix.columns:
        logging.info('Calculating distance for clone {}'.format(cluster_id))
        cell_clone_corr[cluster_id] = cell_cn_matrix.corrwith(
            clone_cn_matrix[cluster_id], method=distance_method)

    distances = pd.DataFrame(cell_clone_corr)
    distances.columns.name = 'cluster_id'
    distances = distances.stack().rename(distance_metric).reset_index()

    return distances


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
        pearsonr_distances.set_index(['cell_id', 'cluster_id']),
        spearmanr_distances.set_index(['cell_id', 'cluster_id']),
        cityblock_distances.set_index(['cell_id', 'cluster_id']),
    ], axis=1).reset_index()

    return clone_cell_distances


def plot_clones(cn_data, cluster_col, plots_prefix):
    plot_data = cn_data.copy()
    bin_filter = (plot_data['gc'] <= 0) | (plot_data['copy'].isnull())
    plot_data.loc[bin_filter, 'state'] = 0
    plot_data.loc[plot_data['copy'] > 5, 'copy'] = 5.
    plot_data.loc[plot_data['copy'] < 0, 'copy'] = 0.

    fig = plt.figure(figsize=(15, 2))
    scgenome.cnplot.plot_cluster_cn_matrix(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    fig.savefig(plots_prefix + '_clone_cn.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'copy', cluster_field_name=cluster_col, raw=True)
    fig.savefig(plots_prefix + '_raw_cn.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    fig.savefig(plots_prefix + '_cn_state.pdf', bbox_inches='tight')


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


def infer_clones(library_ids, sample_ids, results_prefix, local_storage_directory):
    """ Run clonal inference on a set of libraries 
    """
    tantalus_api = dbclients.tantalus.TantalusApi()

    raw_data_filename = results_prefix + '_raw_data.pickle'

    logging.info('retrieving cn data')
    cn_data, metrics_data, image_feature_data = memoize(
        raw_data_filename,
        retrieve_cn_data,
        tantalus_api,
        library_ids,
        sample_ids,
        local_storage_directory,
        results_prefix,
    )

    # TODO: Remove temporary fixup
    if 'total_mapped_reads_hmmcopy' not in metrics_data:
         metrics_data['total_mapped_reads_hmmcopy'] = metrics_data['total_mapped_reads']
    elif metrics_data['total_mapped_reads_hmmcopy'].isnull().any():
        fix_read_count = metrics_data['total_mapped_reads_hmmcopy'].isnull()
        metrics_data.loc[fix_read_count, 'total_mapped_reads_hmmcopy'] = (
            metrics_data.loc[fix_read_count, 'total_mapped_reads'])

    clones_filename = results_prefix + '_clusters.pickle'
    logging.info('calculating clusters')
    shape_check = cn_data.shape
    logging.info('cn_data shape {}'.format(shape_check))
    clusters, filter_metrics = memoize(
        clones_filename,
        calculate_clusters,
        cn_data,
        metrics_data,
        results_prefix,
    )
    assert cn_data.shape == shape_check

    cell_clone_distances_filename = results_prefix + '_cell_cluster_distances.pickle'
    cell_clone_distances = memoize(
        cell_clone_distances_filename,
        calculate_cell_clone_distances,
        cn_data,
        clusters,
        results_prefix,
    )

    final_clusters_filename = results_prefix + '_final_clusters.pickle'
    final_clusters = memoize(
        final_clusters_filename,
        finalize_clusters,
        cn_data,
        metrics_data,
        clusters,
        filter_metrics,
        cell_clone_distances,
        results_prefix,
    )


@click.group()
def infer_clones_cmd():
    pass


@infer_clones_cmd.command('singlelib')
@click.argument('library_id')
@click.argument('sample_id')
@click.argument('results_prefix')
@click.argument('local_storage_directory')
def infer_clones_singlelib_cmd(
        library_id,
        sample_id,
        results_prefix,
        local_storage_directory,
):
    library_ids = [library_id]
    sample_ids = [sample_id]

    infer_clones(
        library_ids,
        sample_ids,
        results_prefix,
        local_storage_directory,
    )


@infer_clones_cmd.command('multilib')
@click.argument('library_ids_filename')
@click.argument('sample_ids_filename')
@click.argument('results_prefix')
@click.argument('local_storage_directory')
def infer_clones_multilib_cmd(
        library_ids_filename,
        sample_ids_filename,
        results_prefix,
        local_storage_directory,
):
    library_ids = [l.strip() for l in open(library_ids_filename).readlines()]
    sample_ids = [l.strip() for l in open(sample_ids_filename).readlines()]

    infer_clones(
        library_ids,
        sample_ids,
        results_prefix,
        local_storage_directory,
    )


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones_cmd()
