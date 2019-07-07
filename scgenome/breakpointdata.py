import seaborn
import logging

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import wgs_analysis.plots.rearrangement


index_cols = ['sample_id', 'library_id', 'prediction_id']

bp_index_cols = [
    'chromosome_1', 'strand_1', 'position_1',
    'chromosome_2', 'strand_2', 'position_2',
]


def load_breakpoint_data(pseudobulk):
    """ Load breakpoints and add annotations.
    """

    breakpoint_data = pseudobulk.load_breakpoint_data()
    breakpoint_count_data = pseudobulk.load_breakpoint_count_data()

    # Calculate cell counts
    cell_counts = (
        breakpoint_count_data
        .query('read_count > 0')
        .drop_duplicates(index_cols + ['cell_id'])
        .groupby(index_cols).size().rename('num_cells')
        .reset_index())

    breakpoint_data = breakpoint_data.merge(cell_counts, how='left')
    assert not breakpoint_data['num_cells'].isnull().any()

    return breakpoint_data, breakpoint_count_data


def filter_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
        num_split_threshold=1,
        template_length_min_threshold=250,
        num_cells_threshold=2,
        num_unique_reads_threshold=2,
    ):
    """ Filter breakpoint data and breakpoint counts
    """

    num_initial_breakpoints = len(breakpoint_data.index)

    filtered_predictions = []

    # Filter by number of split reads
    if num_split_threshold is not None:
        filtered_num_split = breakpoint_data.query(
            'num_split < {}'.format(num_split_threshold))
        logging.info('Filtering {} of {} breakpoints by num_split < {}'.format(
            len(filtered_num_split.index), num_initial_breakpoints, num_split_threshold))
        filtered_predictions.append(filtered_num_split[index_cols].drop_duplicates())

    # Filter by predicted sequence length
    if template_length_min_threshold is not None:
        filtered_template_length_min = breakpoint_data.query(
            'template_length_min < {}'.format(template_length_min_threshold))
        logging.info('Filtering {} of {} breakpoints by template_length_min < {}'.format(
            len(filtered_template_length_min.index), num_initial_breakpoints, template_length_min_threshold))
        filtered_predictions.append(filtered_template_length_min[index_cols].drop_duplicates())

    # Filter breakpoints by num cells in which they are detected
    if num_cells_threshold is not None:
        filtered_cell_counts = breakpoint_data.query(
            'num_cells < {}'.format(num_cells_threshold))
        logging.info('Filtering {} of {} breakpoints by num_cells < {}'.format(
            len(filtered_cell_counts.index), num_initial_breakpoints, num_cells_threshold))
        filtered_predictions.append(filtered_cell_counts[index_cols].drop_duplicates())

    # Filter breakpoints by total read counts
    if num_unique_reads_threshold is not None:
        filtered_num_unique_reads = breakpoint_data.query(
            'num_unique_reads < {}'.format(num_unique_reads_threshold))
        logging.info('Filtering {} of {} breakpoints by num_unique_reads < {}'.format(
            len(filtered_num_unique_reads.index), num_initial_breakpoints, num_unique_reads_threshold))
        filtered_predictions.append(filtered_num_unique_reads[index_cols].drop_duplicates())

    if len(filtered_predictions) > 0:
        filtered_predictions = pd.concat(filtered_predictions, ignore_index=True).drop_duplicates()
        logging.info('Filtering {} of {} breakpoints'.format(
            len(filtered_predictions.index), num_initial_breakpoints))
        breakpoint_data = breakpoint_data[
            ~breakpoint_data.set_index(index_cols).index.isin(filtered_predictions.set_index(index_cols).index)]
        logging.info('Retained {} of {} breakpoints'.format(
            len(breakpoint_data.index), num_initial_breakpoints))

    else:
        logging.info('No filtering applied')

    breakpoint_count_data = breakpoint_count_data.merge(
        breakpoint_data[index_cols + bp_index_cols].drop_duplicates())

    return breakpoint_data, breakpoint_count_data


def plot_library_portrait(breakpoint_data, figures_prefix):
    """ Plot a library level portrait of breakpoint features.
    """

    breakpoint_data['log_num_reads'] = np.log10(breakpoint_data['num_reads'])
    breakpoint_data['log_num_split'] = np.log10(breakpoint_data['num_split'])
    breakpoint_data['log_num_cells'] = np.log10(breakpoint_data['num_cells'])

    # Plot per library feature distributions
    features = [
        ('log_num_reads', 'log discordant read counts'),
        ('log_num_split', 'log split read counts'),
        ('template_length_min', 'prediction sequence length'),
        ('homology', 'sequence homology'),
        ('log_num_cells', 'log number of cells'),
    ]

    order = breakpoint_data['library_id'].unique()
    hue_order = breakpoint_data['rearrangement_type'].unique()

    fig = plt.figure(figsize=(8, 12))
    ax = fig.add_subplot(len(features) + 1, 1, 1)
    plot_data = (
        breakpoint_data
        .groupby(['library_id', 'rearrangement_type'])
        .size().rename('count').reset_index())
    seaborn.barplot(
        ax=ax, data=plot_data,
        x='library_id', y='count', hue='rearrangement_type',
        order=order, hue_order=hue_order)
    ax.set_title(f'Counts by rearrangement type across libraries')
    ax.set_xticklabels([])
    ax.set_xlabel('', visible=False)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    for idx, (col, desc) in enumerate(features):
        ax = fig.add_subplot(len(features) + 1, 1, idx + 2)
        seaborn.boxplot(
            ax=ax, data=breakpoint_data,
            x='library_id', y=col, hue='rearrangement_type',
            order=order, hue_order=hue_order)
        ax.set_title(f'Distribution of {desc} across libraries')
        if idx < len(features) - 1:
            ax.set_xticklabels([])
            ax.set_xlabel('', visible=False)
        ax.legend().set_visible(False)

    plt.tight_layout()
    fig.savefig(figures_prefix + 'library_features.pdf', bbox_inches='tight')

    # Plot rearrangement type distribution across the genome per library
    num_libraries = len(breakpoint_data['library_id'].unique())
    fig = plt.figure(figsize=(10, 3 * num_libraries))
    for idx, (library_id, library_breakpoints) in enumerate(breakpoint_data.groupby('library_id', observed=True)):
        breakends = wgs_analysis.plots.rearrangement.create_breakends(
            library_breakpoints, data_cols=['rearrangement_type'])
        ax = fig.add_subplot(num_libraries, 1, idx + 1)
        wgs_analysis.plots.rearrangement.chromosome_type_plot(
            ax, breakends)
        ax.set_title(f'Chromosome types for library {library_id}')
    plt.tight_layout()
    fig.savefig(figures_prefix + 'chromosome_types.pdf', bbox_inches='tight')

    # Plot adjacent density of breakends across the genome per library
    fig = plt.figure(figsize=(10, 3 * num_libraries))
    for idx, (library_id, library_breakpoints) in enumerate(breakpoint_data.groupby('library_id')):
        breakends = wgs_analysis.plots.rearrangement.create_breakends(
            library_breakpoints)
        breakends['chrom'] = breakends['chromosome']
        breakends['coord'] = breakends['position']
        breakends = wgs_analysis.annotation.position.annotate_adjacent_distance(breakends)
        ax = fig.add_subplot(num_libraries, 1, idx + 1)
        wgs_analysis.plots.snv.snv_adjacent_distance_plot(
            ax, breakends)
        ax.set_title(f'Breakend adjacent distances for library {library_id}')
    plt.tight_layout()
    fig.savefig(figures_prefix + 'adjacent_distance.pdf', bbox_inches='tight')


def plot_breakpoint_clustering(breakpoint_data, breakpoint_count_data, clusters, figures_prefix):
    """ Plot breakpoint cluster figures.
    """
    plot_data = breakpoint_count_data.merge(clusters)
    plot_data = (
        plot_data.groupby(bp_index_cols + ['cluster_id'])[['read_count']]
        .sum().unstack(fill_value=None).stack().reset_index())

    g = seaborn.factorplot(y='read_count', x='cluster_id', kind='box', data=plot_data, color='0.75', size=4)
    g.fig.savefig(figures_prefix + 'cluster_id_read_counts.pdf', bbox_inches='tight')

    plot_data['is_present'] = (plot_data['read_count'] > 0) * 1
    plot_data = plot_data.set_index(bp_index_cols + ['cluster_id'])['is_present'].unstack(fill_value=None)
    mask = plot_data.isnull()
    plot_data = plot_data.fillna(0)

    g = seaborn.clustermap(plot_data, mask=mask, rasterized=True, figsize=(4, 12))
    g.fig.savefig(figures_prefix + 'cluster_map.pdf', bbox_inches='tight')

