import logging
import scipy.stats
import seaborn
import umap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scgenome.cncluster
import scgenome.cnplot
import scgenome.cnfilter


def calculate_umap_hdbscan_clones(cn_data, metrics_data, results_prefix=None):
    """ Cluster copy number data.
    """
    logging.info('creating copy number matrix')
    cn = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0)
    )

    logging.info('clustering copy number')
    clusters = scgenome.cncluster.umap_hdbscan_cluster(cn)

    num_clusters = len(clusters['cluster_id'].unique())

    fig = plt.figure(figsize=(4, 4))
    scgenome.cncluster.plot_umap_clusters(plt.gca(), clusters)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'cn_umap.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(num_clusters/4, 4))
    seaborn.barplot(x='cluster_id', y='count', data=clusters.groupby('cluster_id').size().rename('count').reset_index())
    if results_prefix is not None:
        fig.savefig(results_prefix + 'clone_size.pdf', bbox_inches='tight')

    plot_data = (clusters
        .merge(metrics_data[['cell_id', 'sample_id']].drop_duplicates())
        .groupby(['cluster_id', 'sample_id']).size().unstack(fill_value=0).T)
    fig = plt.figure(figsize=(2, num_clusters/4))
    seaborn.heatmap(plot_data.T, annot=True, fmt="d", linewidths=.5)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'clone_sample_size.pdf', bbox_inches='tight')

    logging.info('merging clusters')
    cn_data = cn_data.merge(clusters[['cell_id', 'cluster_id']].drop_duplicates())

    plot_clones(cn_data, 'cluster_id', plots_prefix=results_prefix)

    return clusters


def calculate_kmeans_clones(cn_data, metrics_data, results_prefix=None, min_k=2, max_k=100, min_cluster_size=None):
    """ Cluster copy number data.
    """
    logging.info('creating copy number matrix')
    cn = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0)
    )

    logging.info('clustering copy number')
    clusters = scgenome.cncluster.kmeans_cluster(cn, min_k=min_k, max_k=max_k)
    clusters = clusters.merge(clusters.groupby('cluster_id').size().rename('cluster_size').reset_index())

    if min_cluster_size is not None:
        clusters = clusters.query(f"cluster_size >= {min_cluster_size}")

    embedding = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        metric='euclidean',
    ).fit_transform(cn.fillna(0).values.T)

    clusters = clusters.merge(pd.DataFrame({
        'cell_id': cn.columns,
        'umap1': embedding[:, 0], 'umap2': embedding[:, 1]
    }))

    num_clusters = len(clusters['cluster_id'].unique())

    fig = plt.figure(figsize=(4, 4))
    scgenome.cncluster.plot_umap_clusters(plt.gca(), clusters)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'cn_umap.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(num_clusters/4, 4))
    seaborn.barplot(x='cluster_id', y='count', data=clusters.groupby('cluster_id').size().rename('count').reset_index())
    if results_prefix is not None:
        fig.savefig(results_prefix + 'clone_size.pdf', bbox_inches='tight')

    plot_data = (clusters
        .merge(metrics_data[['cell_id', 'sample_id']].drop_duplicates())
        .groupby(['cluster_id', 'sample_id']).size().unstack(fill_value=0).T)
    fig = plt.figure(figsize=(2, num_clusters/4))
    seaborn.heatmap(plot_data.T, annot=True, fmt="d", linewidths=.5)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'clone_sample_size.pdf', bbox_inches='tight')

    logging.info('merging clusters')
    cn_data = cn_data.merge(clusters[['cell_id', 'cluster_id']].drop_duplicates())

    plot_clones(cn_data, 'cluster_id', plots_prefix=results_prefix)

    return clusters


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
    fig.savefig(results_prefix + 'breakpoint_hist.pdf', bbox_inches='tight')

    breakpoint_data = breakpoint_data[breakpoint_data['breakpoints'] <= max_breakpoints]

    clusters = clusters[clusters['cluster_id'].isin(breakpoint_data['cluster_id'])]

    return clusters


def finalize_clusters(
        cn_data, metrics_data, clusters, filter_metrics,
        cell_clone_distances, results_prefix,
        is_original_cluster_mean_threshold=0.5,
        cluster_size_threshold=50):
    """ Generate finalized filtered clusters
    """

    # Calculate the cluster assignment based on correlation
    correlation_metric = 'pearsonr'
    correlation_cluster = (
        cell_clone_distances
        .set_index(['cluster_id', 'cell_id'])[correlation_metric]
        .unstack().idxmin().rename(correlation_metric + 'cluster_id').reset_index())

    # Calculate which cells cluster assignment matches highest correlation cluster
    cluster_annotation = cell_clone_distances.merge(clusters[['cell_id', 'cluster_id']])
    cluster_annotation = cluster_annotation.merge(correlation_cluster)
    cluster_annotation['is_original'] = (cluster_annotation['cluster_id'] == cluster_annotation[correlation_metric + 'cluster_id'])

    # Plot the cityblock distance distribution of each cluster separated
    # by whether the assigned cluster equals the correlation cluster
    plot_metric = 'cityblock'
    g = seaborn.factorplot(
        x='cluster_id', y=plot_metric,
        hue='is_original', kind='strip',
        dodge=True, data=cluster_annotation, aspect=3)
    g.fig.savefig(results_prefix + 'cluster_cityblock_distance.pdf', bbox_inches='tight')

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

    # Assign clusters to s phase cells
    if metrics_data['is_s_phase'].any():
        # Assign s phase cells to the cluster they are most correlated with
        cell_filtered_clone_distances = cell_clone_distances.merge(
            cluster_annotation[['cluster_id']].drop_duplicates())
        s_phase_cluster = cell_filtered_clone_distances.merge(
            metrics_data.query('is_s_phase == True')[['cell_id']])
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

    else:
        # No s-phase cells
        final_clusters = cluster_annotation[['cell_id', 'cluster_id']].copy()

    assert not final_clusters['cell_id'].duplicated().any()

    # Plotting
    #
    
    # Plot final clusters heatmap
    logging.info('plotting clusters to {}*'.format(results_prefix + 'final_'))
    plot_cn_data = cn_data.merge(
        final_clusters[['cell_id', 'cluster_id']])
    plot_clones(plot_cn_data, 'cluster_id', results_prefix + 'final_')

    # Plot s phase proportions
    metrics_data['is_s_phase'] = metrics_data['is_s_phase'].astype(bool)
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
    fig.savefig(results_prefix + 'clone_s_phase.pdf', bbox_inches='tight')

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
            .groupby(['chr', 'start', 'end', 'cluster_id'], observed=True)
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


def plot_clones(cn_data, cluster_col, plots_prefix=None):
    plot_data = cn_data.copy()
    bin_filter = (plot_data['gc'] <= 0) | (plot_data['copy'].isnull())
    plot_data.loc[bin_filter, 'state'] = 0
    plot_data.loc[plot_data['copy'] > 5, 'copy'] = 5.
    plot_data.loc[plot_data['copy'] < 0, 'copy'] = 0.

    logging.info('Plotting cluster cn matrix')

    num_clusters = len(cn_data[cluster_col].unique())
    fig = plt.figure(figsize=(15, num_clusters/8))
    scgenome.cnplot.plot_cluster_cn_matrix(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'clone_cn.pdf', bbox_inches='tight')

    logging.info('Plotting cell cn matrix')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'copy', cluster_field_name=cluster_col, raw=True)
    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'raw_cn.pdf', bbox_inches='tight')

    logging.info('Plotting cell cn matrix figure')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'cn_state.pdf', bbox_inches='tight')

    logging.info('Generating clone data')

    clone_cn_data = (
        cn_data
            .groupby(['chr', 'start', 'end', 'cluster_id'], observed=True)
            .agg({'copy': np.mean, 'state': np.median})
            .reset_index()
    )
    clone_cn_data['state'] = clone_cn_data['state'].astype(float).round().astype(int)

    logging.info('Plotting clone data')

    num_clusters = len(clone_cn_data['cluster_id'].unique())
    fig = plt.figure(figsize=(20, 4 * num_clusters))
    idx = 1
    for cluster_id, plot_data in clone_cn_data.groupby('cluster_id'):
        ax = fig.add_subplot(num_clusters, 1, idx)
        scgenome.cnplot.plot_cell_cn_profile(
            ax, plot_data, 'copy', 'state')
        ax.set_ylabel(f'Clone {cluster_id} Total CN')
        idx += 1
    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'total_cn_profiles.pdf', bbox_inches='tight')


def calculate_mitotic_errors(
        cn_data,
        clusters,
        results_prefix,
        chromosome_state_diff_threshold=0.75,
        include_y_chrom=False,
    ):

    # Calculate raw and integer copy state per cluster
    clone_cn_data = (
    cn_data
        .merge(clusters)
        .groupby(['chr', 'start', 'end', 'cluster_id'], observed=True)
        .agg({'copy': np.mean, 'state': np.median})
        .reset_index()
    )
    clone_cn_data['state'] = clone_cn_data['state'].round().astype(int)

    # Create clone / cell cn table 
    clone_cell_cn = cn_data.merge(clusters[['cell_id', 'cluster_id']])
    clone_cell_cn = clone_cell_cn.merge(
        clone_cn_data[['chr', 'start', 'end', 'cluster_id', 'state']].rename(columns={'state': 'clone_cn'}),
        on=['chr', 'start', 'end', 'cluster_id'])
    clone_cell_cn['state_diff'] = clone_cell_cn['state'] - clone_cell_cn['clone_cn']
    clone_cell_cn['bin_size'] = clone_cell_cn['end'] - clone_cell_cn['start'] + 1

    # Calculate proportion of each chromosome for each cell that has
    # a different state from the clone (mean_diff)
    size_state_diff = clone_cell_cn.set_index('state_diff', append=True).groupby(['chr', 'cell_id', 'state_diff'], observed=True)['bin_size'].sum().rename('state_size').reset_index()
    size_total = clone_cell_cn.groupby(['chr', 'cell_id'])['bin_size'].sum().astype(float).rename('total_size').reset_index()
    mean_state_diff = size_state_diff.merge(size_total)
    mean_state_diff['mean_diff'] = mean_state_diff['state_size'] / mean_state_diff['total_size']

    # Filtering and threshold at 0.75 of a chromosome as a mis-segregation
    if not include_y_chrom:
        mean_state_diff = mean_state_diff[mean_state_diff['chr'] != 'Y']
    mean_state_diff = mean_state_diff[mean_state_diff['state_diff'] != 0]
    mean_state_diff = mean_state_diff[mean_state_diff['mean_diff'] > chromosome_state_diff_threshold]

    # Count per state diff
    fig = plt.figure(figsize=(3, 2))
    plot_data = mean_state_diff.groupby('state_diff').size().rename('count').reset_index()
    seaborn.barplot(x='state_diff', y='count', data=plot_data)
    fig.savefig(results_prefix + 'misseg_state_diff_counts.pdf', bbox_inches='tight')

    # Count per state diff
    fig = plt.figure(figsize=(4, 2))
    ax = fig.add_subplot(111)
    plot_data = mean_state_diff.groupby('cell_id').size().rename('chr_count').reset_index()
    plot_data = plot_data.groupby('chr_count').size().rename('cell_count').reset_index()
    plot_data['proportion'] = plot_data['cell_count'] / len(clone_cell_cn['cell_id'].unique())
    plot_data = plot_data.query('chr_count > 0')
    chr_counts = range(1, plot_data['chr_count'].max() + 1)
    seaborn.barplot(ax=ax, x='chr_count', y='proportion', data=plot_data, order=chr_counts, color='0.75')
    ax.set_xlabel('Num. chromosomes')
    ax.set_ylabel('Prop. cells')
    seaborn.despine(trim=True)
    fig.savefig(results_prefix + 'misseg_state_diff_proportions.pdf', bbox_inches='tight')

    # Count per chromosome
    chromosomes = [str(a) for a in range(1, 23)] + ['X']
    if include_y_chrom:
        chromosomes.append('Y')
    fig = plt.figure(figsize=(7, 2.5))
    ax = fig.add_subplot(111)
    plot_data = mean_state_diff.query('state_diff > 0').groupby('chr').size().rename('count').reset_index()
    seaborn.barplot(x='chr', y='count', data=plot_data, order=chromosomes, color=seaborn.desaturate('red', 0.75))
    plot_data = mean_state_diff.query('state_diff < 0').groupby('chr').size().rename('count').reset_index()
    plot_data['count'] = -plot_data['count']
    seaborn.barplot(x='chr', y='count', data=plot_data, order=chromosomes, color=seaborn.desaturate('blue', 0.75))
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Count')
    seaborn.despine(trim=True)
    ax.set_yticklabels([int(abs(a)) for a in ax.get_yticks()])
    fig.savefig(results_prefix + 'misseg_chr_counts.pdf', bbox_inches='tight')

    return mean_state_diff

