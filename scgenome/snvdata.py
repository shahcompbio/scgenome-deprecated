import os
import logging
import seaborn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv


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


