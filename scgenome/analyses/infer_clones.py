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
import scgenome.utils
import scgenome.dataimport
import scgenome.cncluster
import scgenome.cnplot
import scgenome.snvdata
import scgenome.snpdata
import scgenome.breakpointdata
import scgenome.snvphylo

import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv
import wgs_analysis.annotation.position

import dollo
import dollo.tasks

import dbclients.tantalus
from dbclients.basicclient import NotFoundError
import datamanagement.transfer_files


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


# TODO: thresholds
museq_score_threshold = None
strelka_score_threshold = None
snvs_num_cells_threshold = 2
snvs_sum_alt_threshold = 4
is_original_cluster_mean_threshold = 0.5
cluster_size_threshold = 50

# SA1135
# museq_score_threshold = 0.5
# strelka_score_threshold = -np.inf

cn_bin_size = 500000

results_storage_name = 'singlecellblob_results'


def retrieve_cn_data(tantalus_api, library_ids, sample_ids, local_cache_directory, results_prefix):
    hmmcopy_results, hmmcopy_tickets = scgenome.dataimport.search_hmmcopy_analyses(tantalus_api, library_ids)

    results = scgenome.dataimport.import_cn_data(
        hmmcopy_tickets,
        local_cache_directory,
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


def calculate_cluster_allele_cn(cn_data, allele_data, clusters, results_prefix):
    """ Infer allele and cluster specific copy number from haplotype allele counts
    """
    clone_cn_state = (
        cn_data.merge(clusters[['cell_id', 'cluster_id']])
        .groupby(['chr', 'start', 'end', 'cluster_id'])['state']
        .median().astype(int).reset_index())

    clone_cn_copy = (
        cn_data.merge(clusters[['cell_id', 'cluster_id']])
        .groupby(['chr', 'start', 'end', 'cluster_id'])['copy']
        .mean().reset_index())

    clone_cn_data = clone_cn_state.merge(clone_cn_copy)

    clone_cn_data['total_cn'] = clone_cn_data['state']

    allele_data = allele_data.rename(columns={
        'chromosome': 'chr',
        'total': 'total_counts_sum',
        'allele_1': 'allele_1_sum',
        'allele_2': 'allele_2_sum',
    })

    allele_cn = scgenome.snpdata.infer_allele_cn(clone_cn_data, allele_data)

    allele_data['maf'] = (
        np.minimum(allele_data['allele_1_sum'], allele_data['allele_2_sum']) /
        allele_data['total_counts_sum'].astype(float))

    for cluster_id in clusters['cluster_id'].unique():
        cluster_allele_data = allele_data.query('cluster_id == {}'.format(cluster_id))
        cluster_allele_cn = allele_cn.query('cluster_id == {}'.format(cluster_id))

        fig = plt.figure(figsize=(14, 3))
        ax = fig.add_subplot(111)
        scgenome.snpdata.plot_vaf_cn_profile(
            ax, cluster_allele_data, cluster_allele_cn)
        fig.savefig(results_prefix + f'cluster_{cluster_id}_allele_cn.pdf')

    return allele_cn


def infer_allele_cn(allele_data, clone_cn_data, results_prefix):
    """ Infer allele and cluster specific copy number from haplotype allele counts
    """
    import IPython; IPython.embed(); raise

    data = allele_counts.rename(columns={
        'chromosome': 'chr',
        'total': 'total_counts_sum',
        'allele_1': 'allele_1_sum',
        'allele_2': 'allele_2_sum',
    })

    hap_data = data.copy()

    hap_data = hap_data[hap_data['total_counts_sum'] > 5].copy()

    hap_data['maf'] = (
        np.minimum(hap_data['allele_1_sum'], hap_data['allele_2_sum']) /
        hap_data['total_counts_sum'].astype(float))

    allele_cn = scgenome.snpdata.infer_allele_cn(clone_cn_data, data)

    allele_cn['maf'] = (
        np.minimum(allele_cn['allele_1_sum'], allele_cn['allele_2_sum']) /
        allele_cn['total_counts_sum'].astype(float))

    allele_cn.query('chr == "4"').query('total_counts_sum > 100').head()


def retrieve_pseudobulk_data(ticket_id, clusters, local_cache_directory):
    tantalus_api = dbclients.tantalus.TantalusApi()

    results = scgenome.dataimport.search_pseudobulk_results(tantalus_api, ticket_id)
    sample_libraries = scgenome.dataimport.get_pseudobulk_sample_libraries(tantalus_api, ticket_id)

    dataset_filepaths = datamanagement.transfer_files.cache_dataset(
        tantalus_api,
        results['id'],
        'resultsdataset',
        results_storage_name,
        local_cache_directory,
    )

    snv_data = scgenome.snvdata.load_snv_data(dataset_filepaths)

    allele_data = scgenome.snpdata.load_haplotype_allele_counts(dataset_filepaths)
    allele_data = scgenome.snpdata.calculate_cluster_allele_counts(allele_data, clusters)

    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.load_breakpoint_data(
        dataset_filepaths, sample_libraries)

    return snv_data, allele_data, breakpoint_data, breakpoint_count_data


def run_mutation_signature_analysis(snv_data, results_prefix):
    """
    Run a mutation signature analysis
    """
    sigs, sig_prob = wgs_analysis.snvs.mutsig.load_signature_probabilities()
    total_alt_counts = (
        snv_data.groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts']
        .sum().astype(int).rename('total_alt_counts').reset_index())
    snv_data = snv_data.merge(total_alt_counts)

    total_total_counts = (
        snv_data.groupby(['chrom', 'coord', 'ref', 'alt'])['total_counts']
        .sum().astype(int).rename('total_total_counts').reset_index())
    snv_data = snv_data.merge(total_total_counts)

    snv_data = snv_data[snv_data['total_alt_counts'] > 0]
    snv_data['num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] > 1, 'num_cells_class'] = '2-5'
    snv_data.loc[snv_data['num_cells'] > 5, 'num_cells_class'] = '6-20'
    snv_data.loc[snv_data['num_cells'] > 20, 'num_cells_class'] = '>20'

    # Per cell count class signatures
    snv_sig_data = snv_data[snv_data['tri_nucleotide_context'].notnull()]
    snv_sig_data = snv_sig_data.merge(sigs)

    # Simple filter for variant sample presence
    snv_sig_data = snv_sig_data[snv_sig_data['alt_counts'] > 0]

    snv_sig_data2 = snv_sig_data.copy()
    snv_sig_data2['num_cells_class'] = 'All'
    snv_sig_data2 = pd.concat([snv_sig_data, snv_sig_data2])

    # Generate signature distributions for cell count classes
    sample_sig = wgs_analysis.snvs.mutsig.fit_sample_signatures(
        snv_sig_data2.drop_duplicates(
            ['chrom', 'coord', 'ref', 'alt', 'num_cells_class']),
        sig_prob, 'num_cells_class')

    fig = wgs_analysis.snvs.mutsig.plot_signature_heatmap(sample_sig)
    seaborn.despine(trim=True)
    fig.savefig(results_prefix + '_num_cell_class_mutsig.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(4, 4))
    sample_sig.loc['All', :].sort_values().iloc[-10:,].plot(kind='bar')
    seaborn.despine(trim=True)
    fig.savefig(results_prefix + '_top10_mutsig.pdf', bbox_inches='tight')


def run_bulk_snv_analysis(snv_data, results_prefix):
    # Write high impact SNVs to a csv table
    high_impact = (snv_data.query('effect_impact == "HIGH"').query('is_cosmic == True')
        [[
            'chrom', 'coord', 'ref', 'alt',
            'gene_name', 'effect', 'effect_impact',
            'is_cosmic', 'museq_score', 'strelka_score',
        ]]
        .drop_duplicates())
    high_impact = high_impact.merge(
        snv_data.groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_data.groupby(['chrom', 'coord', 'ref', 'alt'])['ref_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_data.groupby(['chrom', 'coord', 'ref', 'alt'])['total_counts'].sum().reset_index())
    high_impact.to_csv(results_prefix + '_snvs_high_impact.csv')

    # Calculate cell counts
    cell_counts = (
        snv_data
        .query('alt_counts > 0')
        .drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
        .groupby(['chrom', 'coord', 'ref', 'alt']).size().rename('num_cells')
        .reset_index())
    fig = plt.figure(figsize=(4, 4))
    cell_counts['num_cells'].hist(bins=50)
    fig.savefig(results_prefix + '_snv_cell_counts.pdf', bbox_inches='tight')
    snv_data = snv_data.merge(cell_counts, how='left')
    assert not snv_data['num_cells'].isnull().any()

    # Calculate total alt counts for each SNV
    sum_alt_counts = (
        snv_data
        .groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts']
        .sum().rename('sum_alt_counts')
        .reset_index())
    fig = plt.figure(figsize=(4, 4))
    sum_alt_counts['sum_alt_counts'].hist(bins=50)
    fig.savefig(results_prefix + '_snv_alt_counts.pdf', bbox_inches='tight')

    # Filter SNVs by num cells in which they are detected
    filtered_cell_counts = cell_counts.query(
        'num_cells >= {}'.format(snvs_num_cells_threshold))
    logging.info('Filtering {} of {} SNVs by num_cells >= {}'.format(
        len(filtered_cell_counts.index), len(cell_counts.index), snvs_num_cells_threshold))
    snv_data = snv_data.merge(filtered_cell_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    # Filter SNVs by total alt counts
    filtered_sum_alt_counts = sum_alt_counts.query(
        'sum_alt_counts >= {}'.format(snvs_sum_alt_threshold))
    logging.info('Filtering {} of {} SNVs by num_cells >= {}'.format(
        len(filtered_sum_alt_counts.index), len(sum_alt_counts.index), snvs_sum_alt_threshold))
    snv_data = snv_data.merge(filtered_sum_alt_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    run_mutation_signature_analysis(snv_data, results_prefix)

    # Plot adjacent distance of SNVs
    fig = plt.figure(figsize=(10, 3))
    snv_data = wgs_analysis.annotation.position.annotate_adjacent_distance(snv_data)
    wgs_analysis.plots.snv.snv_adjacent_distance_plot(plt.gca(), snv_data)
    fig.savefig(results_prefix + '_snv_adjacent_density.pdf', bbox_inches='tight')

    return snv_data


def run_snv_phylogenetics(snv_data, allele_cn, clusters, results_prefix):
    """ Run the SNV phylogenetic analysis.
    """
    snv_log_likelihoods = scgenome.snvphylo.compute_snv_log_likelihoods(
        snv_data, allele_cn, clusters)

    ml_tree, tree_annotations = scgenome.snvphylo.compute_dollo_ml_tree(
        snv_log_likelihoods)

    import IPython; IPython.embed(); raise
    snv_matrix['vaf'] = snv_matrix['alt_counts'] / snv_matrix['total_counts']
    snv_matrix['alt_counts'] = snv_matrix['alt_counts'].clip_upper(10)
    snv_matrix['is_present'] = (snv_matrix['alt_counts'] > 0) * 1
    snv_matrix['is_absent'] = (snv_matrix['alt_counts'] == 0) * 1
    snv_matrix['is_het'] = (snv_matrix['alt_counts'] < 0.99 * snv_matrix['total_counts']) * snv_matrix['is_present']
    snv_matrix['is_hom'] = (snv_matrix['alt_counts'] >= 0.99 * snv_matrix['total_counts']) * snv_matrix['is_present']
    snv_matrix['state'] = snv_matrix['is_hom'] * 3 + snv_matrix['is_het'] * 2 + snv_matrix['is_absent']
    snv_presence_matrix = snv_matrix.set_index(['chrom', 'coord', 'cluster_id'])['is_present'].unstack(fill_value=0)
    print(snv_presence_matrix.shape)

    g = seaborn.clustermap(snv_presence_matrix, rasterized=True, row_cluster=True, figsize=(4, 12))


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


class Memoizer(object):
    def __init__(self, cache_prefix):
        self.cache_prefix = cache_prefix
    def __call__(self, name, func, *args, **kwargs):
        cache_filename = self.cache_prefix + name + '.pickle'
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


def infer_clones(library_ids, sample_ids, pseudobulk_ticket, results_prefix, local_cache_directory):
    """ Run clonal inference on a set of libraries 
    """
    tantalus_api = dbclients.tantalus.TantalusApi()

    memoizer = Memoizer(results_prefix + '_')

    logging.info('retrieving cn data')
    cn_data, metrics_data, image_feature_data = memoizer(
        'raw_data',
        retrieve_cn_data,
        tantalus_api,
        library_ids,
        sample_ids,
        local_cache_directory,
        results_prefix,
    )

    # TODO: Remove temporary fixup
    if 'total_mapped_reads_hmmcopy' not in metrics_data:
         metrics_data['total_mapped_reads_hmmcopy'] = metrics_data['total_mapped_reads']
    elif metrics_data['total_mapped_reads_hmmcopy'].isnull().any():
        fix_read_count = metrics_data['total_mapped_reads_hmmcopy'].isnull()
        metrics_data.loc[fix_read_count, 'total_mapped_reads_hmmcopy'] = (
            metrics_data.loc[fix_read_count, 'total_mapped_reads'])

    logging.info('calculating clusters')
    shape_check = cn_data.shape
    logging.info('cn_data shape {}'.format(shape_check))
    clusters, filter_metrics = memoizer(
        'clusters',
        calculate_clusters,
        cn_data,
        metrics_data,
        results_prefix,
    )
    assert cn_data.shape == shape_check

    cell_clone_distances = memoizer(
        'cell_cluster_distances',
        calculate_cell_clone_distances,
        cn_data,
        clusters,
        results_prefix,
    )

    final_clusters = memoizer(
        'final_clusters',
        finalize_clusters,
        cn_data,
        metrics_data,
        clusters,
        filter_metrics,
        cell_clone_distances,
        results_prefix,
    )

    snv_data, allele_data, breakpoint_data, breakpoint_count_data = memoizer(
        'pseudobulk_data',
        retrieve_pseudobulk_data,
        pseudobulk_ticket,
        final_clusters,
        local_cache_directory,
    )

    allele_cn = memoizer(
        'allele_cn',
        calculate_cluster_allele_cn,
        cn_data,
        allele_data,
        clusters,
        results_prefix,
    )
    
    filtered_snv_data = memoizer(
        'snv_filtering',
        run_bulk_snv_analysis,
        snv_data,
        results_prefix,
    )

    snv_phylogeny = memoizer(
        'snv_phylogeny',
        run_snv_phylogenetics,
        filtered_snv_data,
        allele_cn,
        final_clusters,
        results_prefix, 
    )


@click.group()
def infer_clones_cmd():
    pass


@infer_clones_cmd.command('singlelib')
@click.argument('library_id')
@click.argument('sample_id')
@click.argument('pseudobulk_ticket')
@click.argument('results_prefix')
@click.argument('local_cache_directory')
def infer_clones_singlelib_cmd(
        library_id,
        sample_id,
        pseudobulk_ticket,
        results_prefix,
        local_cache_directory,
):
    library_ids = [library_id]
    sample_ids = [sample_id]

    infer_clones(
        library_ids,
        sample_ids,
        pseudobulk_ticket,
        results_prefix,
        local_cache_directory,
    )


@infer_clones_cmd.command('multilib')
@click.argument('library_ids_filename')
@click.argument('sample_ids_filename')
@click.argument('pseudobulk_ticket')
@click.argument('results_prefix')
@click.argument('local_cache_directory')
def infer_clones_multilib_cmd(
        library_ids_filename,
        sample_ids_filename,
        pseudobulk_ticket,
        results_prefix,
        local_cache_directory,
):
    library_ids = [l.strip() for l in open(library_ids_filename).readlines()]
    sample_ids = [l.strip() for l in open(sample_ids_filename).readlines()]

    infer_clones(
        library_ids,
        sample_ids,
        pseudobulk_ticket,
        results_prefix,
        local_cache_directory,
    )


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones_cmd()
