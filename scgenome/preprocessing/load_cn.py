import anndata as ad

import scgenome.loaders.qc


def read_dlp_hmmcopy(alignment_results_dir, hmmcopy_results_dir, annotation_results_dir, sample_ids=None, additional_hmmcopy_reads_cols=None):
    """ Read hmmcopy results from the DLP pipeline.

    Parameters
    ------
    alignment_results_dir (str):
        dlp pipeline alignment results directory
    hmmcopy_results_dir (str):
        dlp pipeline hmmcopy results directory
    annotation_results_dir (str):
        dlp pipeline annotation results directory
    sample_ids (list):
        sample ids to load
    additional_hmmcopy_reads_cols (list):
        per bin metrics to load

    Returns
    ------
    AnnData
        An instantiated AnnData Object.
    """

    results = scgenome.loaders.qc.load_qc_results(
        alignment_results_dir,
        hmmcopy_results_dir,
        annotation_results_dir,
        sample_ids=sample_ids,
        additional_hmmcopy_reads_cols=additional_hmmcopy_reads_cols,
    )

    metrics_data = results['annotation_metrics']
    cn_data = results['hmmcopy_reads']

    cn_matrix = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])[['reads', 'copy', 'state']]
            .unstack(level='cell_id')
            .transpose())

    bin_data = (
        cn_data[['chr', 'start', 'end', 'gc', 'map']]
            .drop_duplicates()
            .set_index(['chr', 'start', 'end'])
            .reindex(cn_matrix.loc['reads'].columns))

    cell_data = (
        metrics_data
            .set_index(['cell_id'])
            .reindex(cn_matrix.loc['reads'].index))

    adata = ad.AnnData(
        cn_matrix.loc['reads'],
        obs=cell_data,
        var=bin_data,
        layers={
            'copy': cn_matrix.loc['copy'],
            'state': cn_matrix.loc['state'],
        },
    )

    return adata
