import logging
import seaborn

import numpy as np
import pandas as pd
import scipy.stats
import scipy.misc
import scipy.special
import matplotlib.pyplot as plt
import matplotlib

import dollo.tasks
import dollo.run
import wgs_analysis.plots.trees

import scgenome.snvphylo


def annotate_copy_number(pos, seg, columns=['major', 'minor'], sample_col='sample_id'):
    """ Annotate positions with segment specific data 
    """
    results = []

    for sample_id in seg[sample_col].unique():
        sample_pos = pos[pos[sample_col] == sample_id]
        sample_seg = seg[seg[sample_col] == sample_id]

        for chrom in seg['chr'].unique():
            _pos = sample_pos[sample_pos['chrom'] == chrom]
            _seg = sample_seg[sample_seg['chr'] == chrom]

            results.append(find_overlapping_segments(_pos, _seg, columns))

    return pd.concat(results)


def find_overlapping_segments(pos, seg, columns):
    """ Find positions that are contained within segments
    """
    seg = seg.sort_values(['start', 'end'])

    if seg.duplicated(['start', 'end']).any():
        raise ValueError('duplicate columns')

    start_idx = np.searchsorted(seg['start'].values, pos['coord'].values) - 1
    end_idx = np.searchsorted(seg['end'].values, pos['coord'].values)

    mask = (start_idx == end_idx)

    results = pos.copy()

    for col in columns:
        results[col] = np.nan
        results.loc[mask, col] = seg[col].iloc[end_idx[mask]].values

    return results


def snv_hierarchical_clustering_figure(snv_data, clusters):
    """ Simple hierarhical clustering figure for SNVs
    """
    snv_matrix = snv_data.merge(clusters)

    snv_matrix = (
        snv_matrix.groupby(
            ['chrom', 'coord', 'ref', 'alt', 'cluster_id'],
            as_index=True, observed=True)[['alt_counts', 'ref_counts']]
        .sum().unstack().fillna(0).astype(int).stack().reset_index())
    snv_matrix['total_counts'] = snv_matrix['ref_counts'] + snv_matrix['alt_counts']

    snv_matrix['vaf'] = snv_matrix['alt_counts'] / snv_matrix['total_counts']
    snv_matrix['alt_counts'] = snv_matrix['alt_counts'].clip(upper=10)
    snv_matrix['is_present'] = (snv_matrix['alt_counts'] > 0) * 1
    snv_matrix['is_absent'] = (snv_matrix['alt_counts'] == 0) * 1
    snv_matrix['is_het'] = (snv_matrix['alt_counts'] < 0.99 * snv_matrix['total_counts']) * snv_matrix['is_present']
    snv_matrix['is_hom'] = (snv_matrix['alt_counts'] >= 0.99 * snv_matrix['total_counts']) * snv_matrix['is_present']
    snv_matrix['state'] = snv_matrix['is_hom'] * 3 + snv_matrix['is_het'] * 2 + snv_matrix['is_absent']
    snv_presence_matrix = snv_matrix.set_index(['chrom', 'coord', 'cluster_id'])['is_present'].unstack(fill_value=0)

    logging.info(f'snv matrix with shape {snv_presence_matrix.shape}, memory {snv_presence_matrix.memory_usage().sum()}')

    # KLUDGE: currently recursion in dendrograms
    # breaks with large datasets
    import sys
    sys.setrecursionlimit(10000)

    g = seaborn.clustermap(snv_presence_matrix, rasterized=True, row_cluster=True, figsize=(4, 12))

    return g.fig


def compute_snv_log_likelihoods(snv_data, allele_cn, clusters):
    """ Compute log likelihoods of presence absence for SNVs
    """
    snv_matrix = snv_data.merge(clusters)

    snv_matrix = (
        snv_matrix.groupby(
            ['chrom', 'coord', 'ref', 'alt', 'cluster_id'],
            as_index=True, observed=True)[['alt_counts', 'ref_counts']]
        .sum().unstack().fillna(0).astype(int).stack().reset_index())
    snv_matrix['total_counts'] = snv_matrix['ref_counts'] + snv_matrix['alt_counts']

    snv_matrix['variant_id'] = snv_matrix.apply(
        lambda row: ':'.join(row[['chrom', 'coord', 'ref', 'alt']].astype(str).values),
        axis=1).astype('category')

    # TODO: this should be moved
    allele_cn['total_cn'] = allele_cn['total_cn'].astype(int)
    allele_cn['minor_cn'] = allele_cn['minor_cn'].astype(int)
    allele_cn['major_cn'] = allele_cn['major_cn'].astype(int)
    allele_cn = allele_cn[[
        'chr', 'start', 'end', 'cluster_id',
        'total_cn', 'minor_cn', 'major_cn',
    ]].drop_duplicates()

    # Merge segment copy number into SNV table
    snv_log_likelihoods = annotate_copy_number(
        snv_matrix, allele_cn,
        columns=['major_cn', 'minor_cn', 'total_cn'],
        sample_col='cluster_id')

    snv_log_likelihoods = compute_log_likelihoods(
        snv_log_likelihoods)

    return snv_log_likelihoods


def compute_log_likelihoods(df, error_rate=1e-3):
    """ Compute the presence absence log likelihood of an SNV
    """
    df['log_likelihood_absent'] = df.apply(calculate_likelihood_absent, axis=1, args=(error_rate,))
    df['log_likelihood_present'] = df.apply(calculate_likelihood_present, axis=1, args=(error_rate,))

    return df


def calculate_likelihood_absent(row, e_s):
    return log_likelihood_absent(
        e_s,
        row['alt_counts'],
        row['alt_counts'] + row['ref_counts'],
    )


def calculate_likelihood_present(row, e_s):
    return log_likelihood_present(
        row['major_cn'],
        row['major_cn'] + row['minor_cn'],
        e_s,
        row['alt_counts'],
        row['ref_counts'] + row['alt_counts'],
    )

def log_binomial_pdf(x, n, p):
    return scipy.stats.binom.logpmf(x, n, p)


def log_likelihood_absent(e_s, n_v, n_t):
    return log_binomial_pdf(n_v, n_t, e_s)


def log_likelihood_present(c_m, c_t, e_s, n_v, n_t):
    if c_m == 0:
        return log_likelihood_absent(e_s, n_v, n_t)

    conditional_log_likelihoods = []

    for c_v in np.arange(1., c_m + 1., 1.):
        allele_ratio = c_v / c_t
        r = (1 - e_s) * allele_ratio + e_s * (1 - allele_ratio)
        conditional_log_likelihoods.append(log_binomial_pdf(n_v, n_t, r))

    return scipy.special.logsumexp(conditional_log_likelihoods)


def compute_dollo_ml_tree(snv_log_likelihoods, leaf_name_groups=None):
    """ Compute the ML tree under the dollo model of SNV evolution
    """
    trees = dollo.tasks.create_trees(
        snv_log_likelihoods,
        sample_col='cluster_id',
        leaf_name_groups=leaf_name_groups,
    )

    results_table = dollo.tasks.compute_tree_log_likelihoods_mp(
        snv_log_likelihoods, trees,
        sample_col='cluster_id', variant_col='variant_id')

    ml_tree_id = results_table.set_index('tree_id')['log_likelihood'].idxmax()
    tree = trees[ml_tree_id]
    loss_prob = results_table.set_index('tree_id').loc[ml_tree_id, 'loss_prob']

    tree_annotations = dollo.run.annotate_posteriors(
        snv_log_likelihoods, tree, loss_prob=loss_prob,
        sample_col='cluster_id', variant_col='variant_id')

    return tree, tree_annotations


def plot_dollo_ml_tree(tree, nodes):
    """ Plot a tree with branchlengths as snv origin counts
    """
    leaf_order = []
    for leaf in tree.leaves:
        leaf.plot_id = leaf.name
        leaf_order.append(leaf.name)

    origin_counts = nodes.groupby('node')['ml_origin'].sum()

    for node in tree.nodes:
        node.origin_count = origin_counts[node.label]

    loss_counts = nodes.groupby('node')['ml_loss'].sum()

    width = 1 + 0.5 * float(len(list(tree.leaves)))
    fig = plt.figure(figsize=(width/1.5, 6))

    ax = fig.add_subplot(111)

    def func(x, pos):
        s = '{:0,d}'.format(int(x))
        return s
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(func))

    wgs_analysis.plots.trees.plot_tree(ax, tree, landscape=False, flip=True, branch_length_attr='origin_count', leaf_name_attr='plot_id')

    ax.set_ylabel('SNV count')

    plt.tight_layout()


def run_snv_phylogenetics(snv_count_data, allele_cn, clusters, results_prefix):
    """ Run the SNV phylogenetic analysis.
    """
    snv_log_likelihoods = scgenome.snvphylo.compute_snv_log_likelihoods(
        snv_count_data, allele_cn, clusters)

    ml_tree, tree_annotations = scgenome.snvphylo.compute_dollo_ml_tree(
        snv_log_likelihoods)

    return ml_tree, tree_annotations

