import scipy
import numpy as np
import pandas as pd


def filter_cells(cn, scores, threshold=0.75):
    cn.set_index(['chr', 'start', 'end', 'width'], inplace=True)

    scores = scores[
        (scores['joint'] >= threshold) &
        (scores['cell_call'].isin(['C1', 'C2']))
    ]

    for rm_cond in ['gDNA', 'GM', 'NCC', 'NTC']:
        mask = ~scores['experimental_condition'].str.contains(rm_cond)
        scores = scores[mask]

    good_cells = scores.loc[
        scores['cell_id'].isin(cn.columns.tolist()), 'cell_id'
    ].tolist()
    good_cells.sort()

    return cn[good_cells].reset_index()


def filter_qdnaseq_bins(cn, blacklist):
    blacklist.rename(columns={'chromosome': 'chr'}, inplace=True)

    cn = cn.merge(blacklist, how='left')

    # residual cutoff is 4 * madDiff(residual)
    cn = cn[
        (
            (
                ~pd.isnull(cn['residual']) &
                (cn['residual'].abs() <= 0.0696924)
            ) |
            pd.isnull(cn['residual'])
        ) &
        (cn['blacklist'] == 0)
    ]

    rm_cols = ['bases', 'gc', 'mappability', 'blacklist', 'residual', 'use']
    for rm_col in rm_cols:
        del cn[rm_col]

    return cn


def remove_contiguous_duplicate_bins(cn):
    cn.sort_values(by=['chr', 'start', 'end'], inplace=True)
    cn.set_index(['start', 'end', 'width'], inplace=True)
    cn = cn[(cn.shift() != cn).any(axis=1)]
    cn = cn.reset_index().set_index(['chr', 'start', 'end', 'width'])
    cn.reset_index(inplace=True)
    return cn


def calc_prop_hom_del(states):
    cndist = states.value_counts()
    cndist = cndist / cndist.sum()
    if 0 not in cndist:
        return 0
    return cndist[0]


def calculate_filter_metrics(
        metrics_data,
        cn_data,
        quality_score_threshold=0.75,
        read_count_threshold=500000,
        prop_hom_del_pval_threshold=0.01,
        copy_state_diff_threshold=1.,
    ):
    """ Calculate additional filtering values and add to metrics.
    """
    metrics_data['filter_quality'] = (metrics_data['quality'] > quality_score_threshold)
    metrics_data['filter_reads'] = (metrics_data['total_mapped_reads_hmmcopy'] > read_count_threshold)

    # Calculate proportion homozygous deletion state
    prop_hom_del = cn_data.groupby('cell_id')['state'].apply(calc_prop_hom_del).rename('prop_hom_del').reset_index()
    metrics_data = metrics_data.merge(prop_hom_del, how='left')
    metrics_data['prop_hom_del'] = metrics_data['prop_hom_del'].fillna(0)

    # Calculate p value for each proportion assuming beta fit
    a, b, loc, scale = scipy.stats.beta.fit(metrics_data['prop_hom_del'].values)
    metrics_data['prop_hom_del_pval'] = 1.0 - scipy.stats.beta.cdf(
        metrics_data['prop_hom_del'].values, a, b, loc, scale)
    metrics_data['filter_prop_hom_del'] = (metrics_data['prop_hom_del_pval'] > prop_hom_del_pval_threshold)

    # Calculate separation between predicted and normalized copy number
    copy_state_diff = cn_data[['cell_id', 'copy', 'state']].copy()
    copy_state_diff['copy_state_diff'] = np.absolute(copy_state_diff['copy'] - copy_state_diff['state'])
    copy_state_diff = (copy_state_diff[['cell_id', 'copy_state_diff']]
        .dropna().groupby('cell_id')['copy_state_diff']
        .mean().reset_index().dropna())
    metrics_data = metrics_data.merge(copy_state_diff)
    metrics_data['filter_copy_state_diff'] = (metrics_data['copy_state_diff'] < copy_state_diff_threshold)

    # Filter s phase column
    metrics_data['filter_is_s_phase'] = ~(metrics_data['is_s_phase'].fillna(False))

    return metrics_data
