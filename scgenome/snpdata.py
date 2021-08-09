import matplotlib
import seaborn
import logging
import scipy.special
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from hmmlearn._hmmc import _viterbi
from scipy.stats import binom
from matplotlib import collections  as mc

import scgenome.refgenome
import scgenome.utils
import scgenome.cnplot


# TODO: possibly deprecated
def get_clone_snp_count(het_counts, cluster_df):
    cell_ids = het_counts['cell_ids']
    cell_ids.set_index('cell_id', inplace=True)

    snp_data = None

    cluster_df = cluster_df.merge(cell_ids.reset_index()[['cell_id']].drop_duplicates())

    print(cluster_df.shape[0])
    for cell_id, cluster_id in cluster_df[['cell_id', 'cluster_id']].values:
        print('.', end=' ')
        cell_index = cell_ids.loc[cell_id, 'cell_index']
        df = het_counts.select('/allele_counts', where='cell_index=={}'.format(cell_index))
        df['cluster_id'] = cluster_id
        df['chromosome'] = df['chromosome'].astype(str)
        if snp_data is None:
            snp_data = df
        else:
            snp_data = pd.concat([snp_data, df], ignore_index=True)
            snp_data = snp_data.groupby(['chromosome', 'coord', 'cluster_id'], observed=True)[['ref_counts', 'alt_counts']].sum().reset_index()

    snp_data['total_counts'] = snp_data['ref_counts'] + snp_data['alt_counts']
    snp_data['vaf'] = snp_data['alt_counts'] / snp_data['total_counts'].astype(float)

    return snp_data


# TODO: possibly deprecated
def calculate_haplotype_allele_counts(snp_data, haps, bin_size):
    haps['swap'] = (haps['allele'] == haps['allele_id'])
    haps = haps[['chromosome', 'position', 'hap_label', 'swap']].drop_duplicates().rename(columns={'position': 'coord'})

    data = snp_data.merge(haps)
    data['allele_1'] = np.where(data['swap'], data['ref_counts'], data['alt_counts'])
    data['allele_2'] = np.where(~data['swap'], data['ref_counts'], data['alt_counts'])
    data['bin'] = (data['coord'] / bin_size).astype(int)

    # Aggregated by bin and haplotype block
    data2 = (
        data.groupby(['cluster_id', 'chromosome', 'hap_label', 'bin'], observed=True)[['allele_1', 'allele_2', 'total_counts', 'coord']]
        .aggregate({'allele_1': sum, 'allele_2': sum, 'total_counts': sum, 'coord': (np.min, np.max, len)}).reset_index())
    data2.columns = ['_'.join(col).strip().rstrip('_') for col in data2.columns.values]
    data2['maf'] = np.minimum(data2['allele_1_sum'], data2['allele_2_sum']) / data2['total_counts_sum'].astype(float)
    data2['start'] = data2['bin'] * bin_size + 1
    data2['end'] = (data2['bin'] + 1) * bin_size

    hap_data = data2[[
        'chromosome', 'start', 'end', 'hap_label',
        'cluster_id', 'allele_1_sum', 'allele_2_sum', 'total_counts_sum'
    ]].rename(columns={'chromosome': 'chr'})

    return hap_data


def plot_vaf_cn_profile(ax, hap_data, allele_cn):
    """ Plot genome wide VAF and predicted CN state for haplotype alleles.
    """
    allele_cn = allele_cn.query('total_cn != 0')

    chrom_info = plot_vaf_profile(
        ax, hap_data, 'maf',
        size_field_name='total_counts_sum',
        size_scale=100.,
    )

    allele_cn['cn_ratio'] = allele_cn['minor_cn'] / allele_cn['total_cn']

    allele_cn = allele_cn[allele_cn['chr'].isin(scgenome.refgenome.info.chromosomes)]

    allele_cn.set_index('chr', inplace=True)
    allele_cn['chromosome_start'] = scgenome.refgenome.info.chromosome_start
    allele_cn.reset_index(inplace=True)

    allele_cn['start'] = allele_cn['start'] + allele_cn['chromosome_start']
    allele_cn['end'] = allele_cn['end'] + allele_cn['chromosome_start']

    lines = np.zeros((allele_cn.shape[0], 2, 2))
    lines[:, 0, 0] = allele_cn['start'].values
    lines[:, 1, 0] = allele_cn['end'].values
    lines[:, 0, 1] = allele_cn['cn_ratio'].values
    lines[:, 1, 1] = allele_cn['cn_ratio'].values

    lc = mc.LineCollection(lines, colors='b', linewidths=1.)
    ax.add_collection(lc)

    return chrom_info


def plot_vaf_profile(ax, cn_data, value_field_name, cn_field_name=None, size_field_name=None, size_scale=1.):
    """ Plot genome wide VAF profile for haplotype alleles
    """
    plot_data = cn_data.copy()
    plot_data = plot_data[plot_data['chr'].isin(scgenome.refgenome.info.chromosomes)]

    plot_data.set_index('chr', inplace=True)
    plot_data['chromosome_start'] = scgenome.refgenome.info.chromosome_start
    plot_data.reset_index(inplace=True)

    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    c = '0.75'
    cmap = None
    if cn_field_name is not None:
        c = plot_data[cn_field_name]
        cmap = scgenome.cnplot.get_cn_cmap(plot_data[cn_field_name].values)

    s = 1
    if size_field_name is not None:
        s = plot_data[size_field_name] / size_scale

    ax.scatter(
        plot_data['start'], plot_data[value_field_name], c=c, s=s, alpha=0.1, cmap=cmap)

    ax.set_xlim((-0.5, scgenome.refgenome.info.chromosome_end.max()))
    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(scgenome.refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.set_yticks([0., 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(scgenome.refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(scgenome.refgenome.info.chromosomes))

    seaborn.despine(ax=ax, offset=10, trim=True)

    chrom_info = (
        plot_data[['chr', 'chromosome_start']]
        .drop_duplicates()
        .set_index('chr')['chromosome_start'])

    return chrom_info


def infer_allele_cn(clone_cn_data, hap_data, loh_error_rate=0.01):
    """ HMM inference of clone and allele specific copy number based on haplotype
    allele read counts.
    """
    cn = clone_cn_data.merge(
        hap_data,
        on=['chr', 'start', 'end', 'cluster_id'],
        how='left',
    ).fillna(0).rename(columns={'integer_copy_number': 'total_cn'})

    cn['total_cn'] = cn['total_cn'].astype(int).astype(float)

    total_cn = cn['total_cn'].values

    minor_cn = np.tile(np.arange(0, 10, 1), (cn.shape[0], 1))
    n = np.tile(cn['total_counts_sum'].astype(int), (minor_cn.shape[1], 1)).T
    x = np.tile(cn[['allele_1_sum', 'allele_2_sum']].min(axis=1).astype(int), (minor_cn.shape[1], 1)).T
    p = minor_cn / total_cn[:, np.newaxis]

    # Add an error term for the loh state
    # Rational: we can think of each states probability of a minor allele
    # read as being beta distributed around the expected calculated as
    # the ratio of minor to total copy number.  For most states the expectation
    # will be equal to the copy number ratio except for the loh state which will
    # have the configurable offset from 0 governed by the loh_error_rate term
    p[:, 0] = loh_error_rate

    l = binom.logpmf(x, n, p)
    l = np.logaddexp(l, np.log(0.01))

    # Set the likelihood to be very low for impossible states
    l[minor_cn > total_cn[:, np.newaxis]] = -1000.

    # Set the likelihood to a low value that is still greater than the
    # impossible state value if the likelihood is nan.  This will primarily
    # catch the case where minor copy number equals total copy number, which
    # is in general not a valid solution unless total copy number is 0
    l[np.isnan(l)] = -100.

    minor_cn_data = (pd.DataFrame(l, index=cn.set_index(['cluster_id', 'chr', 'start', 'end', 'hap_label']).index)
        .groupby(level=[0, 1, 2, 3]).sum())

    minor_cn = pd.DataFrame(index=minor_cn_data.index)
    minor_cn['minor_cn'] = None

    for clone_id in cn['cluster_id'].unique():
        framelogprob = minor_cn_data.loc[clone_id].values

        N_chain = framelogprob.shape[0]
        N_states = framelogprob.shape[1]

        seq, prob = _viterbi(
            N_chain,
            N_states,
            np.log(np.ones(N_states) / N_states),
            np.log(np.eye(N_states) * 1e4 + np.ones((N_states, N_states))),
            framelogprob)

        minor_cn.loc[clone_id, 'minor_cn'] = seq

    minor_cn = minor_cn.reset_index()
    minor_cn = minor_cn.merge(
        cn[['chr', 'start', 'end', 'cluster_id', 'total_cn']].drop_duplicates(),
        how='left')
    minor_cn['total_cn'] = minor_cn['total_cn'].astype(int)
    minor_cn['major_cn'] = minor_cn['total_cn'] - minor_cn['minor_cn']

    assert minor_cn['minor_cn'].notnull().any()
    assert minor_cn['major_cn'].notnull().any()
    assert minor_cn['total_cn'].notnull().any()

    return minor_cn


def calculate_cluster_allele_counts(allele_data, clusters, cn_bin_size):
    """ Calculate allele specific haplotype allele counts per cluster
    """
    index_cols = [
        'chromosome',
        'start',
        'end',
        'hap_label',
    ]

    # Create allele 1 and 2 counts matrix
    allele_data = allele_data.set_index(index_cols + ['cell_id', 'allele_id'])['readcount'].unstack(fill_value=0)
    assert '0' in allele_data.columns and '1' in allele_data.columns
    allele_data.rename(columns={'0': 'allele_1', '1': 'allele_2'}, inplace=True)
    allele_data.reset_index(inplace=True)

    # Merge clusters and redo categoricals
    scgenome.utils.union_categories(
        [allele_data, clusters],
        cat_cols=['cell_id'])
    allele_data = allele_data.merge(clusters[['cell_id', 'cluster_id']])
    allele_data = allele_data.groupby(
        index_cols + ['cluster_id'],
        as_index=True, observed=True)[['allele_1', 'allele_2']].sum().reset_index()

    allele_data['total'] = allele_data['allele_1'] + allele_data['allele_2']
    allele_data['start'] = (allele_data['start'] / cn_bin_size).astype(int) * cn_bin_size + 1
    allele_data['end'] = allele_data['start'] + cn_bin_size - 1

    return allele_data


def calculate_cluster_allele_cn(
        cn_data, allele_data, clusters,
        total_allele_counts_threshold=6,
        plots_prefix=None,
    ):
    """ Infer allele and cluster specific copy number from haplotype allele counts
    """
    clone_cn_state = (
        cn_data.merge(clusters[['cell_id', 'cluster_id']])
        .groupby(['chr', 'start', 'end', 'cluster_id'], observed=True)['state']
        .median().astype(int).reset_index())

    clone_cn_copy = (
        cn_data.merge(clusters[['cell_id', 'cluster_id']])
        .groupby(['chr', 'start', 'end', 'cluster_id'], observed=True)['copy']
        .mean().reset_index())

    clone_cn_data = clone_cn_state.merge(clone_cn_copy)

    clone_cn_data['total_cn'] = clone_cn_data['state']

    allele_data = allele_data.rename(columns={
        'chromosome': 'chr',
        'total': 'total_counts_sum',
        'allele_1': 'allele_1_sum',
        'allele_2': 'allele_2_sum',
    })

    allele_data = allele_data[allele_data['total_counts_sum'] >= total_allele_counts_threshold]

    allele_cn = scgenome.snpdata.infer_allele_cn(clone_cn_data, allele_data)

    allele_data['maf'] = (
        np.minimum(allele_data['allele_1_sum'], allele_data['allele_2_sum']) /
        allele_data['total_counts_sum'].astype(float))

    num_clusters = len(clusters['cluster_id'].unique())
    fig = plt.figure(figsize=(20, 4 * num_clusters))
    idx = 1
    for cluster_id in clusters['cluster_id'].unique():
        cluster_allele_data = allele_data.query('cluster_id == {}'.format(cluster_id))
        cluster_allele_cn = allele_cn.query('cluster_id == {}'.format(cluster_id))

        ax = fig.add_subplot(num_clusters, 1, idx)
        scgenome.snpdata.plot_vaf_cn_profile(
            ax, cluster_allele_data, cluster_allele_cn)
        ax.set_ylabel(f'Clone {cluster_id} BAF')

        idx += 1

    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'allele_cn_profiles.png', bbox_inches='tight')

    return allele_cn

