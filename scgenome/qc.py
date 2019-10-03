import scipy.stats
import numpy as np
import pandas as pd


def calc_prop_hom_del(states):
    cndist = states.value_counts()
    cndist = cndist / cndist.sum()
    if 0 not in cndist:
        return 0
    return cndist[0]


def qc_cn(metrics_data, cn_data, quality=0.75, n_reads=500000):
    # TODO handle case where total_mapped_reads is name and not total_mapped_reads_hmmcopy
    if quality is not None:
        metrics_data['filter_quality'] = (metrics_data['quality'] > quality)
    if n_reads is not None:
        metrics_data['filter_reads'] = (metrics_data['total_mapped_reads_hmmcopy'] > n_reads)
    print(f"metrics_data.shape {metrics_data.shape}")

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
    # Remove cells with divergence between copy state and norm copy number
    # Remove cells with outlier proportion of homozygous deletion
    bool_filter = (
        #(~metrics_data['is_s_phase']) &
        metrics_data['filter_copy_state_diff'] &
        metrics_data['filter_prop_hom_del']
    )
    if 'is_s_phase' in metrics_data.columns:
      bool_filter = bool_filter & (~metrics_data['is_s_phase'])
    if quality is not None:
        bool_filter = bool_filter & metrics_data['filter_quality']
    if n_reads is not None:
        bool_filter = bool_filter & metrics_data['filter_reads']
    filtered_cells = metrics_data.loc[bool_filter, ['cell_id']]
    #filtered_cells = metrics_data.loc[
    #    (~metrics_data['is_s_phase']) &
    #    metrics_data['filter_quality'] &
    #    metrics_data['filter_reads'] &
    #    metrics_data['filter_copy_state_diff'] &
    #    metrics_data['filter_prop_hom_del'],
    #    ['cell_id']]

    print(f"filtered_cells.shape: {filtered_cells.shape}")

    cn_data = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())

    cn = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0)
    )

    return cn, cn_data


def filt_reads(metrics_data, cn_data, n_reads=None, quality=None):
    metrics_data['filter_quality'] = (metrics_data['quality'] > quality)
    metrics_data['filter_reads'] = (
            metrics_data['total_mapped_reads_hmmcopy'] > n_reads)

    if quality is not None and n_reads is not None:
        bool_filter = (metrics_data['filter_reads'] &
                       metrics_data['filter_quality'])
    elif n_reads is not None:
        bool_filter = metrics_data['filter_reads']
    elif quality is not None:
        bool_filter = metrics_data['filter_quality']
    else:
        raise ValueError("Must provide 1 of n_reads, quality")
    filtered_cells = metrics_data.loc[bool_filter, ['cell_id']]

    cn_data = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())
    return cn_data


def limit_reads(metrics_data, cn_data, n_reads=None, quality=None):
    metrics_data['filter_quality'] = (metrics_data['quality'] > quality)
    metrics_data['filter_reads'] = (
            metrics_data['total_mapped_reads_hmmcopy'] <= n_reads)

    if quality is not None and n_reads is not None:
        bool_filter = (metrics_data['filter_reads'] &
                       metrics_data['filter_quality'])
    elif n_reads is not None:
        bool_filter = metrics_data['filter_reads']
    elif quality is not None:
        bool_filter = metrics_data['filter_quality']
    else:
        raise ValueError("Must provide 1 of n_reads, quality")
    filtered_cells = metrics_data.loc[bool_filter, ['cell_id']]

    cn_data = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())
    return cn_data


