import matplotlib
import seaborn
import scipy.special
import pandas as pd
import numpy as np
from hmmlearn._hmmc import _viterbi
from scipy.stats import binom

import refgenome


def get_clone_snp_count(het_counts, cluster_df):
    cell_ids = het_counts['cell_ids']
    cell_ids.set_index('cell_id', inplace=True)

    snp_data = None

    cluster_df = cluster_df.merge(cell_ids.reset_index()[['cell_id']].drop_duplicates())

    print cluster_df.shape[0]
    for cell_id, cluster_id in cluster_df[['cell_id', 'cluster_id']].values:
        print '.',
        cell_index = cell_ids.loc[cell_id, 'cell_index']
        df = het_counts.select('/allele_counts', where='cell_index=={}'.format(cell_index))
        df['cluster_id'] = cluster_id
        df['chromosome'] = df['chromosome'].astype(str)
        if snp_data is None:
            snp_data = df
        else:
            snp_data = pd.concat([snp_data, df], ignore_index=True)
            snp_data = snp_data.groupby(['chromosome', 'coord', 'cluster_id'])[['ref_counts', 'alt_counts']].sum().reset_index()

    snp_data['total_counts'] = snp_data['ref_counts'] + snp_data['alt_counts']
    snp_data['vaf'] = snp_data['alt_counts'] / snp_data['total_counts'].astype(float)

    return snp_data


def calculate_haplotype_allele_counts(snp_data, haps, bin_size):
    haps['swap'] = (haps['allele'] == haps['allele_id'])
    haps = haps[['chromosome', 'position', 'hap_label', 'swap']].drop_duplicates().rename(columns={'position': 'coord'})

    data = snp_data.merge(haps)
    data['allele_1'] = np.where(data['swap'], data['ref_counts'], data['alt_counts'])
    data['allele_2'] = np.where(~data['swap'], data['ref_counts'], data['alt_counts'])
    data['bin'] = (data['coord'] / bin_size).astype(int)

    # Aggregated by bin and haplotype block
    data2 = (
        data.groupby(['cluster_id', 'chromosome', 'hap_label', 'bin'])[['allele_1', 'allele_2', 'total_counts', 'coord']]
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


def plot_cell_vaf_profile(ax, cn_data, value_field_name, cn_field_name=None, size_field_name=None, size_scale=1.):
    plot_data = cn_data.copy()
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]

    plot_data.set_index('chr', inplace=True)
    plot_data['chromosome_start'] = refgenome.info.chromosome_start
    plot_data.reset_index(inplace=True)

    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    c = '0.75'
    cmap = None
    if cn_field_name is not None:
        c = plot_data[cn_field_name]
        cmap = get_cn_cmap(plot_data[cn_field_name].values)

    s = 1
    if size_field_name is not None:
        s = plot_data[size_field_name] / size_scale

    ax.scatter(
        plot_data['start'], plot_data[value_field_name], c=c, s=s, alpha=0.1, cmap=cmap)

    ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.chromosomes))

    seaborn.despine(offset=10, trim=True)


def infer_allele_cn(clone_cn_data, hap_data):
    cn = clone_cn_data.merge(
        hap_data,
        on=['chr', 'start', 'end', 'cluster_id'],
        how='left',
    ).fillna(0).rename(columns={'integer_copy_number': 'total_cn'})

    cn['total_cn'] = cn['total_cn'].astype(int).astype(float)

    minor_cn = np.tile(np.arange(0, 10, 1), (cn.shape[0], 1))
    n = np.tile(cn['total_counts_sum'].astype(int), (minor_cn.shape[1], 1)).T
    x = np.tile(cn[['allele_1_sum', 'allele_2_sum']].min(axis=1).astype(int), (minor_cn.shape[1], 1)).T
    p = minor_cn / cn['total_cn'].values[:, np.newaxis]
    p += 0.01

    l = binom.logpmf(x, n, p)
    l = np.logaddexp(l, np.log(0.01))
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

    assert minor_cn['minor_cn'].notnull().any()

    allele_cn = cn.merge(minor_cn.reset_index())
    allele_cn['major_cn'] = allele_cn['total_cn'] - allele_cn['minor_cn']

    return allele_cn

