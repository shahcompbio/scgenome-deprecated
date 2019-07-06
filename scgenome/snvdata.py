import os
import logging
import seaborn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import scgenome.utils
import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv


def get_highest_snpeff_effect(snpeff_data):
    """ Select the highest ranked effect from snpeff data.
    """
    ordered_effect_impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

    ordered_effect_impacts = pd.DataFrame({
            'effect_impact': ordered_effect_impacts,
            'effect_impact_rank': range(len(ordered_effect_impacts))})
    ordered_effect_impacts['effect_impact'] = (
        ordered_effect_impacts['effect_impact'].astype(snpeff_data['effect_impact'].dtype))

    snpeff_data = snpeff_data.merge(ordered_effect_impacts)

    index_cols = ['chrom', 'coord', 'ref', 'alt']
    value_cols = ['gene_name', 'effect', 'effect_impact', 'amino_acid_change']

    snpeff_data = (
        snpeff_data[index_cols + value_cols + ['effect_impact_rank']]
        .sort_values(index_cols + ['effect_impact_rank'], ascending=True)
        .groupby(index_cols, sort=False, observed=True)
        .nth(0)
        .reset_index())

    snpeff_data = snpeff_data[[
        'chrom', 'coord', 'ref', 'alt',
        'gene_name', 'effect', 'effect_impact', 'amino_acid_change',
    ]]

    return snpeff_data


categorical_columns = [
    'chrom',
    'ref',
    'alt',
    'gene_name',
    'effect',
    'effect_impact',
    'amino_acid_change',
    'tri_nucleotide_context',
]

default_museq_filter = 0.9
default_strelka_filter = 20.

def get_snv_results(pseudobulk, museq_filter=None, strelka_filter=None):
    """ Collate snv results into a single table.
    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    logging.info('starting load')
    mappability = pseudobulk.load_snv_annotation_data('mappability')
    mappability = mappability[mappability['mappability'] > 0.99]
    mappability['chrom'] = mappability['chrom'].astype(str)

    strelka_results = pseudobulk.load_snv_annotation_data('strelka').rename(columns={'score': 'strelka_score'})
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = pseudobulk.load_snv_annotation_data('museq').rename(columns={'score': 'museq_score'})
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = pseudobulk.load_snv_annotation_data('cosmic_status')
    logging.info('cosmic table with shape {}'.format(cosmic.shape))
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = pseudobulk.load_snv_annotation_data('snpeff')
    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info('snpeff table with shape {}'.format(snpeff.shape))

    tnc = pseudobulk.load_snv_annotation_data('trinuc')

    data = pseudobulk.load_snv_annotation_data('allele_counts')

    logging.info('summing snv counts')
    data = (
        data
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)[['alt_counts', 'ref_counts']]
        .sum().rename(columns={'alt_counts': 'alt_counts_sum', 'ref_counts': 'ref_counts_sum'}).reset_index())
    logging.info('total snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(mappability)
    logging.info('post mappability with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(cosmic, how='left')
    logging.info('post cosmic with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(snpeff, how='left')
    logging.info('post snpeff with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(tnc, how='left')

    data = data.merge(strelka_results, how='left')
    logging.info('post strelka with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(museq_results, how='left')
    logging.info('post museq with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    if museq_filter != -np.inf:
        data = data[data['museq_score'] > museq_filter]
        logging.info('post museq filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    if strelka_filter != -np.inf:
        data = data[data['strelka_score'] > strelka_filter]
        logging.info('post strelka filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    logging.info('finishing load with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    for column in categorical_columns:
        data[column] = data[column].astype('category')

    return data


def load_snv_data(
        pseudobulk,
        museq_filter=None,
        strelka_filter=None,
        num_cells_threshold=None,
        sum_alt_threshold=None,
        figures_prefix=None,
):
    """ Load filtered SNV annotation and count data
    
    Args:
        pseudobulk (PseudobulkData): pseudbulk data to load from
        museq_score_threshold (float, optional): mutationseq score threshold. Defaults to None.
        strelka_score_threshold (float, optional): strelka score threshold. Defaults to None.
        num_cells_threshold (int, optional): minimum number of cells threshold. Defaults to None.
        sum_alt_threshold (int, optional): minimum total alt count threshold. Defaults to None.
        figures_prefix (str, optional): filename prefix for figures. Defaults to None.
    
    Returns:
        pandas.DataFrame, pandas.DataFrame: SNV annotations, SNV counts
    """
    snv_data = get_snv_results(
        pseudobulk,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter)

    assert not snv_data['coord'].isnull().any()

    positions = snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

    snv_count_data = pseudobulk.load_snv_count_data(positions)
    snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
    snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

    assert not snv_count_data['alt_counts'].isnull().any()

    # Calculate cell counts
    cell_counts = (
        snv_count_data
        .query('alt_counts > 0')
        .drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
        .groupby(['chrom', 'coord', 'ref', 'alt']).size().rename('num_cells')
        .reset_index())

    if figures_prefix is not None:
        fig = plt.figure(figsize=(4, 4))
        cell_counts['num_cells'].hist(bins=50)
        fig.savefig(figures_prefix + 'snv_cell_counts.pdf', bbox_inches='tight')

    snv_data = snv_data.merge(cell_counts, how='left')
    assert not snv_data['num_cells'].isnull().any()

    # Calculate total alt counts for each SNV
    sum_alt_counts = (
        snv_count_data
        .groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts']
        .sum().rename('sum_alt_counts')
        .reset_index())

    if figures_prefix is not None:
        fig = plt.figure(figsize=(4, 4))
        sum_alt_counts['sum_alt_counts'].hist(bins=50)
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


def plot_mutation_signatures(snv_data, snv_class_col, results_prefix):
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
    fig.savefig(results_prefix + f'_{snv_class_col}_mutsig.pdf', bbox_inches='tight')

    fig = plt.figure(figsize=(4, 4))
    sample_sig.loc['All', :].sort_values().iloc[-10:,].plot(kind='bar')
    seaborn.despine(trim=True)
    fig.savefig(results_prefix + f'_{snv_class_col}_top10_mutsig.pdf', bbox_inches='tight')


def run_mutation_signature_analysis(snv_data, results_prefix):
    """
    Run a mutation signature analysis
    """

    # Per cell count class signatures
    snv_data = snv_data[snv_data['num_cells'] > 0]
    snv_data['num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] > 1, 'num_cells_class'] = '2-5'
    snv_data.loc[snv_data['num_cells'] > 5, 'num_cells_class'] = '6-20'
    snv_data.loc[snv_data['num_cells'] > 20, 'num_cells_class'] = '>20'

    plot_mutation_signatures(snv_data, 'num_cells_class', results_prefix)

    # Adjacent distance class signatures
    snv_data['adjacent_distance_class'] = 'standard'
    snv_data.loc[snv_data['adjacent_distance'] <= 10000, 'adjacent_distance_class'] = 'hypermutated'

    plot_mutation_signatures(snv_data, 'adjacent_distance_class', results_prefix)


def run_bulk_snv_analysis(snv_data, snv_count_data, filtered_cell_ids, results_prefix):
    # Filter cells
    snv_count_data = snv_count_data.merge(filtered_cell_ids)
    total_alt_counts = snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts'].sum().reset_index()
    snv_data = snv_data.merge(
        total_alt_counts.query('alt_counts > 0')[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    # Write high impact SNVs to a csv table
    high_impact = (snv_data.query('effect_impact == "HIGH"')
        [[
            'chrom', 'coord', 'ref', 'alt',
            'gene_name', 'effect', 'effect_impact',
            'is_cosmic', 'museq_score', 'strelka_score',
        ]]
        .drop_duplicates())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'])['ref_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'])['total_counts'].sum().reset_index())
    high_impact.to_csv(results_prefix + 'snvs_high_impact.csv')

    # Annotate adjacent distance
    snv_data = wgs_analysis.annotation.position.annotate_adjacent_distance(snv_data)

    # Plot adjacent distance of SNVs
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_adjacent_distance_plot(plt.gca(), snv_data)
    fig.savefig(results_prefix + 'snv_adjacent_distance.pdf', bbox_inches='tight')

    # Plot snv count as a histogram across the genome
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_count_plot(plt.gca(), snv_data)
    fig.savefig(results_prefix + 'snv_genome_count.pdf', bbox_inches='tight')

    # Run mutation signature analysis, requires adjacent distance
    run_mutation_signature_analysis(snv_data, results_prefix)


