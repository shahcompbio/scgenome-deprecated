import logging
import numpy as np
import scipy.stats

from anndata import AnnData



def _calc_prop_hom_del(states):
    cndist, cndist_values = np.unique(states, return_counts=True)
    cndist_values = cndist_values / cndist_values.sum()
    if (0 not in cndist) or (0.0 not in cndist):
        return float(0)
    return float(cndist_values[0])


_default_filters = (
    'filter_quality',
    'filter_reads',
    'filter_copy_state_diff',
    'filter_is_s_phase',
    'filter_prop_hom_del'
)

def calculate_filter_metrics(
        adata: AnnData,
        quality_score_threshold=0.75,
        read_count_threshold=500000,
        prop_hom_del_pval_threshold=0.01,
        copy_state_diff_threshold=1.,
        inplace = False,
    ) -> AnnData:
    """ Calculate additional filtering metrics to be used by other filtering methods.

    Parameters
    ----------
    adata : AnnData
        copy number data on which to calculate filter metrics
    quality_score_threshold : float, optional
        The minimum quality to set to keep, by default 0.75
    read_count_threshold : int, optional
        The minimum total mapped reads from hmmcopy to set for keeping, by default 500000
    prop_hom_del_pval_threshold : float, optional
        Minimum p-value of the proportion of homozygous deletion state to keep, by default 0.01
    copy_state_diff_threshold : [type], optional
        Minimum copy-state difference threshold to set to keep, by default 1.
    inplace : bool, optional
        Whether to modify passed in AnnData, by default False

    Returns
    -------
    AnnData
        AnnData with modified obs if not inplace, otherise, None
    
    Properties Changed
    ------
    AnnData.obs.filter_quality
    AnnData.obs.filter_reads
    AnnData.obs.filter_prop_hom_del
    AnnData.obs.filter_copy_state_diff
    
    If is_s_phase is a property of AnnData
        AnnData.obs.filter_is_s_phase
        
    AnnData.obs.prop_hom_del
    AnnData.obs.prop_hom_del_pval
    AnnData.obs.copy_state_diff
    AnnData.obs.copy_state_diff_mean
    """
    if not inplace:
        ndad = adata.copy()
        
        return calculate_filter_metrics(
            ndad,
            quality_score_threshold,
            read_count_threshold,
            prop_hom_del_pval_threshold,
            copy_state_diff_threshold,
            inplace = True,
        )

    # Filter Quality and Filter Reads
    if 'quality' in adata.obs.columns:
        adata.obs['filter_quality'] = (adata.obs['quality'] > quality_score_threshold)
    else:
        logging.warning("quality is not in AnnData.obs. Skipping filter_quality")
    
    if 'total_mapped_reads_hmmcopy' in adata.obs.columns:
        adata.obs['filter_reads'] = (adata.obs['total_mapped_reads_hmmcopy'] > read_count_threshold)
    else:
        logging.warning("total_mapped_reads_hmmcopy is not in AnnData.obs. Skipping total_mapped_reads_hmmcopy")
    
    # Calculate homozygous deletion state
    adata.obs['prop_hom_del'] = np.apply_along_axis(_calc_prop_hom_del, 1, adata.layers['state'])

    a, b, loc, scale = scipy.stats.beta.fit(adata.obs['prop_hom_del'])
    adata.obs['prop_hom_del_pval'] = 1.0 - scipy.stats.beta.cdf(
        adata.obs['prop_hom_del'], a, b, loc, scale)
    adata.obs['filter_prop_hom_del'] = (adata.obs['prop_hom_del_pval'] > prop_hom_del_pval_threshold)
    adata.obs['prop_hom_del'].fillna(0.0)
    
    # Copy State Difference Filter
    adata.obsm['copy_state_diff'] = np.absolute(adata.layers['copy'] - adata.layers['state'])
    adata.obsm['copy_state_diff_mean'] = adata.obsm['copy_state_diff'].mean(axis=1)

    adata.obs['filter_copy_state_diff'] = (adata.obsm['copy_state_diff_mean'] < copy_state_diff_threshold)

    # Filter s phase column
    if 'is_s_phase' in adata.obs.columns:
        adata.obs['filter_is_s_phase'] = ~(adata.obs['is_s_phase'].fillna(False))
    else:
        logging.warning("No is_s_phase in AnnData.obs. Skipping filter_is_s_phase")

    return adata


def filter_cells(
        adata: AnnData,
        filters = _default_filters,
        inplace = False,
    ) -> AnnData:
    """
    Filter poor quality cells based on the filters provided.

    Parameters
    -------
    adata : AnnData
        AnnData to preform operation with
    filters : list, optional
        Filters to apply. Keeps cells where filters are true, by default _default_filters
    inplace
        Whether to modify passed in AnnData. If False, returns new AnnData.

    Returns
    -------
    AnnData
        filtered copy number data
    """

    if not inplace:
        adata = adata.copy()
        
        return filter_cells(
            adata,
            filters,
            inplace = True,
        )

    # Ensure cnfilter.calculate_filter_metrics has been called
    for filter_option in filters:
        if filter_option not in adata.obs.columns:
            logging.warning(
                f"WARNING: {filter_option} is not found! ",
                "Skipping. Are you sure `scgenome.pp.calculate_filter_metrics` has been called?"
            )
            continue

        adata = adata[adata.obs[filter_option]]

    return adata


