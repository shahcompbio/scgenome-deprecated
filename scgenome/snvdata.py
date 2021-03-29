import os
import logging
import seaborn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv
import wgs_analysis.annotation.position


def filter_snv_data(
        snv_data,
        snv_count_data,
        num_cells_threshold=None,
        sum_alt_threshold=None,
        figures_prefix=None,
    ):
    """ Load filtered SNV annotation and count data

    Args:
        snv_data (pandas.DataFrame): snv annotation data table
        snv_count_data (pandas.DataFrame): snv count data table
        num_cells_threshold (int, optional): minimum number of cells threshold. Defaults to None.
        sum_alt_threshold (int, optional): minimum total alt count threshold. Defaults to None.
        figures_prefix (str, optional): filename prefix for figures. Defaults to None.

    Returns:
        pandas.DataFrame, pandas.DataFrame: SNV annotations, SNV counts
    """
    logging.info('Filtering and annotating SNVs')

    # Calculate cell counts
    cell_counts = (
        snv_count_data
        .query('alt_counts > 0')
        .drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True).size().rename('num_cells')
        .reset_index())

    fig = plt.figure(figsize=(10, 4))
    cell_counts['num_cells'].astype(float).hist(bins=50)
    plt.xlabel('Num. Cells')
    plt.ylabel('Log Num SNVs')
    plt.gca().set_yscale('log')
    if figures_prefix is not None:
        fig.savefig(figures_prefix + 'snv_cell_counts.pdf', bbox_inches='tight')

    snv_data = snv_data.merge(cell_counts, how='left')
    if snv_data['num_cells'].isnull().any():
        num_no_count_snvs = snv_data['num_cells'].isnull().sum()
        logging.warning(f'found {num_no_count_snvs} with no counts in genotyping results')
        snv_data['num_cells'] = snv_data['num_cells'].fillna(0).astype(int)

    # Calculate total alt counts for each SNV
    sum_alt_counts = (
        snv_count_data
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts']
        .sum().rename('sum_alt_counts')
        .reset_index())

    fig = plt.figure(figsize=(10, 4))
    sum_alt_counts['sum_alt_counts'].astype(float).hist(bins=50)
    plt.xlabel('Total Alt. Counts')
    plt.ylabel('Log Num SNVs')
    plt.gca().set_yscale('log')
    if figures_prefix is not None:
        fig.savefig(figures_prefix + 'snv_alt_counts.pdf', bbox_inches='tight')

    # Filter SNVs by num cells in which they are detected
    if num_cells_threshold is not None:
        filtered_cell_counts = cell_counts.query(
            'num_cells >= {}'.format(num_cells_threshold))
        logging.info('Filtering {} of {} SNVs by num_cells >= {}'.format(
            len(filtered_cell_counts.index), len(cell_counts.index), num_cells_threshold))
        snv_data = snv_data.merge(filtered_cell_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    # Filter SNVs by total alt counts
    if sum_alt_threshold is not None:
        filtered_sum_alt_counts = sum_alt_counts.query(
            'sum_alt_counts >= {}'.format(sum_alt_threshold))
        logging.info('Filtering {} of {} SNVs by sum_alt_counts >= {}'.format(
            len(filtered_sum_alt_counts.index), len(sum_alt_counts.index), sum_alt_threshold))
        snv_data = snv_data.merge(filtered_sum_alt_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    snv_count_data = snv_count_data.merge(snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    return snv_data, snv_count_data


def plot_mutation_signatures(snv_data, snv_class_col, figures_prefix=None):
    """
    Infer mutation signature probabilities and plot.
    """

    sigs, sig_prob = wgs_analysis.snvs.mutsig.load_signature_probabilities()

    snv_data = snv_data[snv_data['tri_nucleotide_context'].notnull()]
    snv_data = snv_data.merge(sigs)

    # Add in an all snvs category
    snv_data_all = snv_data.copy()
    snv_data_all[snv_class_col] = 'All'
    snv_data = pd.concat([snv_data, snv_data_all])

    # Generate signature distributions for cell count classes
    sample_sig = wgs_analysis.snvs.mutsig.fit_sample_signatures(
        snv_data.drop_duplicates(
            ['chrom', 'coord', 'ref', 'alt', snv_class_col]),
        sig_prob, snv_class_col)

    fig = wgs_analysis.snvs.mutsig.plot_signature_heatmap(sample_sig)
    seaborn.despine(trim=True)
    if figures_prefix is not None:
        fig.savefig(figures_prefix + f'_{snv_class_col}_mutsig.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(4, 4))
    if len(sample_sig.index) > 0:
        sample_sig.loc['All', :].sort_values().iloc[-10:,].plot(kind='bar')
    seaborn.despine(trim=True)
    if figures_prefix is not None:
        fig.savefig(figures_prefix + f'_{snv_class_col}_top10_mutsig.pdf', bbox_inches='tight')


def run_mutation_signature_analysis(snv_data, figures_prefix=None):
    """
    Run a mutation signature analysis
    """

    # Per cell count class signatures
    snv_data = snv_data[snv_data['num_cells'] > 0]
    snv_data['num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] > 1, 'num_cells_class'] = '2-5'
    snv_data.loc[snv_data['num_cells'] > 5, 'num_cells_class'] = '6-20'
    snv_data.loc[snv_data['num_cells'] > 20, 'num_cells_class'] = '>20'

    plot_mutation_signatures(snv_data, 'num_cells_class', figures_prefix=figures_prefix)

    # Adjacent distance class signatures
    snv_data['adjacent_distance_class'] = 'standard'
    snv_data.loc[snv_data['adjacent_distance'] <= 10000, 'adjacent_distance_class'] = 'hypermutated'

    plot_mutation_signatures(snv_data, 'adjacent_distance_class', figures_prefix=figures_prefix)


def run_bulk_snv_analysis(snv_data, snv_count_data, filtered_cell_ids, results_prefix=None):
    # Filter cells
    snv_count_data = snv_count_data.merge(filtered_cell_ids)
    total_alt_counts = snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts'].sum().reset_index()
    snv_data = snv_data.merge(
        total_alt_counts.query('alt_counts > 0')[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    # Write high impact SNVs to a csv table
    high_impact = snv_data.query('effect_impact == "HIGH"')
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['ref_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['total_counts'].sum().reset_index())
    if results_prefix is not None:
        high_impact.to_csv(results_prefix + 'snvs_high_impact.csv')

    # Annotate adjacent distance
    snv_data = wgs_analysis.annotation.position.annotate_adjacent_distance(snv_data)

    # Plot adjacent distance of SNVs
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_adjacent_distance_plot(plt.gca(), snv_data)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'snv_adjacent_distance.pdf', bbox_inches='tight')

    # Plot snv count as a histogram across the genome
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_count_plot(plt.gca(), snv_data)
    if results_prefix is not None:
        fig.savefig(results_prefix + 'snv_genome_count.pdf', bbox_inches='tight')

    # Run mutation signature analysis, requires adjacent distance
    run_mutation_signature_analysis(snv_data, results_prefix)
