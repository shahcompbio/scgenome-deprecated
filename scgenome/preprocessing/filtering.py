
def _calc_prop_hom_del(states):
    cndist, cndist_values = np.unique(states, return_counts=True)
    cndist_values = cndist_values / cndist_values.sum()
    if (0 not in cndist) or (0.0 not in cndist):
        return float(0)
    return float(cndist_values[0])

default_filters = [
    'filter_quality',
    'filter_reads',
    'filter_copy_state_diff',
    'filter_is_s_phase',
    'filter_prop_hom_del'
]

def calculate_filter_metrics(
        dad,
        quality_score_threshold=0.75,
        read_count_threshold=500000,
        prop_hom_del_pval_threshold=0.01,
        copy_state_diff_threshold=1.,
        inplace = False,
        *args, **kwargs
    ):
    """
    Calculate additional filtering metrics to be used by other filtering methods.
    
    Args
    -------
    dad (DLPAnnData)
        DLPAnnData to preform operation with
    inplace
        Whether to modify passed in DLPAnnData. If False, returns new DLPAnnData.
        
    Optional
    ------
    quality_score_threshold
        The minimum quality to set to keep.
    read_count_threshold
        The minimum total mapped reads from hmmcopy to set for keeping.
    prop_hom_del_pval_threshold
        Minimum p-value of the proportion of homozygous deletion state to keep.
    copy_state_diff_threshold
        Minimum copy-state difference threshold to set to keep.
    
    Returns
    ------
    DLPAnnData with modified obs if not inplace.
    Otherise, None
    
    Properties Changed
    ------
    DLPAnnData.obs.filter_quality
    DLPAnnData.obs.filter_reads
    DLPAnnData.obs.filter_prop_hom_del
    DLPAnnData.obs.filter_copy_state_diff
    
    If is_s_phase is a property of DLPAnnData
        DLPAnnData.obs.filter_is_s_phase
        
    DLPAnnData.obs.prop_hom_del
    DLPAnnData.obs.prop_hom_del_pval
    DLPAnnData.obs.copy_state_diff
    DLPAnnData.obs.copy_state_diff_mean
    """
    if not inplace:
        ndad = dad.copy()
        
        return calculate_filter_metrics(
            ndad,
            quality_score_threshold,
            read_count_threshold,
            prop_hom_del_pval_threshold,
            copy_state_diff_threshold,
            inplace = True,
            *args, **kwargs
        )

    # Filter Quality and Filter Reads
    if 'quality' in dad.obs.columns:
        dad.obs['filter_quality'] = (dad.obs['quality'] > quality_score_threshold)
    else:
        print("WARNING: quality is not in DLPAnnData.obs. Skipping filter_quality")
    
    if 'total_mapped_reads_hmmcopy' in dad.obs.columns:
        dad.obs['filter_reads'] = (dad.obs['total_mapped_reads_hmmcopy'] > read_count_threshold)
    else:
        print("WARNING: total_mapped_reads_hmmcopy is not in DLPAnnData.obs. Skipping total_mapped_reads_hmmcopy")
    
    # Calculate homozygous deletion state
    dad.obs['prop_hom_del'] = np.apply_along_axis(_calc_prop_hom_del, 1, dad.set_X('state').X)

    a, b, loc, scale = scipy.stats.beta.fit(dad.obs['prop_hom_del'])
    dad.obs['prop_hom_del_pval'] = 1.0 - scipy.stats.beta.cdf(
        dad.obs['prop_hom_del'], a, b, loc, scale
    )
    dad.obs['filter_prop_hom_del'] = (dad.obs['prop_hom_del_pval'] > prop_hom_del_pval_threshold)
    dad.obs['prop_hom_del'].fillna(0.0)
    
    # Copy State Difference Filter
    dad.obsm['copy_state_diff'] = np.absolute(dad.set_X('copy').X - dad.set_X('state').X)
    dad.obsm['copy_state_diff_mean'] = dad.obsm['copy_state_diff'].mean(axis=1)

    dad.obs['filter_copy_state_diff'] = (dad.obsm['copy_state_diff_mean'] < copy_state_diff_threshold)

    # Filter s phase column
    if 'is_s_phase' in dad.obs.columns:
        dad.obs['filter_is_s_phase'] = ~(dad.obs['is_s_phase'].fillna(False))
    else:
        print("WARNING: No is_s_phase in DLPAnnData.obs. Skipping filter_is_s_phase")
    

def filter_cells(dad,
                 filters = default_filters,
                 inplace = False,
                 *args, **kwargs):
    """
    Filter poor quality cells based on the filters provided.

    Args
    -------
    dad (DLPAnnData)
        DLPAnnData to preform operation with
    filters [default = default_filters]
        Filters to apply. Keeps cells where filters are true
    inplace
        Whether to modify passed in DLPAnnData. If False, returns new DLPAnnData.
    """

    # Ensure cnfilter.calculate_filter_metrics has been called
    if not inplace:
        dad = dad.copy()
        
        return filter_cells(
                 dad,
                 filters,
                 inplace = True,
                 *args, **kwargs
        )

    for filter_option in filters:
        if filter_option not in fdad.obs.columns:
            print(
                f"WARNING: {filter_option} is not found! ",
                "Skipping. Are you sure `scgenome.pp.calculate_filter_metrics` has been called?"
            )
            continue

        dad = dad[ dad.obs[filter_option] ]




