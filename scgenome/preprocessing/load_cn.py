import logging
import pandas as pd
import anndata as ad
import pyranges as pr

import scgenome.loaders.qc
import csverve

from anndata import AnnData
from pyranges import PyRanges
from typing import Dict
from pandas import DataFrame


def read_dlp_hmmcopy(alignment_results_dir, hmmcopy_results_dir, annotation_results_dir, sample_ids=None, additional_hmmcopy_reads_cols=None) -> AnnData:
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

    return convert_dlp_hmmcopy(metrics_data, cn_data)


def read_dlp_hmmcopy2(reads_filename, metrics_filename, sample_ids=None) -> AnnData:
    """ Read hmmcopy results from the DLP pipeline.

    Parameters
    ------
    reads_filename (str):
        dlp pipeline reads filename
    metrics_filename (str):
        dlp pipeline metrics filename
    sample_ids (list):
        sample ids to load

    Returns
    ------
    AnnData
        An instantiated AnnData Object.
    """

    cn_data = csverve.read_csv(
        reads_filename,
        dtype={
            'cell_id': 'category',
            'sample_id': 'category',
            'library_id': 'category',
            'chr': 'category',
        })

    metrics_data = csverve.read_csv(
        metrics_filename,
        dtype={
            'cell_id': 'category',
            'sample_id': 'category',
            'library_id': 'category',
        })

    scgenome.utils.union_categories([cn_data, metrics_data])

    return convert_dlp_hmmcopy(metrics_data, cn_data)


def convert_dlp_hmmcopy(metrics_data: DataFrame, cn_data: DataFrame) -> AnnData:
    """ Convert hmmcopy pandas dataframes to anndata

    Parameters
    ----------
    metrics_data : DataFrame
        hmmcopy metrics
    cn_data : DataFrame
        hmmcopy reads data

    Returns
    -------
    AnnData
        An instantiated AnnData Object.
    """

    cn_data['bin'] = cn_data['chr'].astype(str) + ':' + cn_data['start'].astype(str) + '-' + cn_data['end'].astype(str)

    cn_matrix = (
        cn_data
            .set_index(['bin', 'cell_id'])[['reads', 'copy', 'state']]
            .unstack(level='cell_id')
            .transpose())

    bin_data = (
        cn_data
            .drop(['cell_id', 'sample_id', 'library_id', 'reads', 'copy', 'state'], axis=1)
            .drop_duplicates(subset=['bin'])
            .set_index(['bin'])
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


def _convert_pyranges(r):
    return r.df.rename(columns={
        'Chromosome': 'chr',
        'Start': 'start',
        'End': 'end',
    })


def _add_bin_index(df):
    df['bin'] = (
        df['chr'].astype(str) + ':' +
        df['start'].astype(str) + '-' +
        df['end'].astype(str))
    return df.set_index('bin')
    

def read_bam_bin_counts(bins: PyRanges, bams: Dict[str, str], excluded: PyRanges = None, **kwargs) -> AnnData:
    """ Count reads in bins from bams

    Parameters
    ----------
    bins : pyranges.PyRanges
        bins in which to count reads
    bams : Dict[Str]
        bam filenames with cell ids as keys
    excluded: PyRanges
        excluded genomic regions to filter reads

    Returns
    -------
    ad.AnnData
        binned read counts
    """

    bin_data = _convert_pyranges(bins)
    bin_data = _add_bin_index(bin_data)

    cn_matrix = {}

    for cell_id, cell_bam in bams.items():
        logging.info(f"reading {cell_bam}")
        bam_data = pr.read_bam(cell_bam, **kwargs)

        if excluded is not None:
            logging.info("excluding reads")
            bam_data = bam_data.intersect(excluded, invert=True)

        logging.info(f"count overlaps")
        bam_data = bam_data.intersect(bins, how='containment')
        read_counts = bins.count_overlaps(bam_data, overlap_col='reads')

        read_counts = _convert_pyranges(read_counts)
        read_counts = _add_bin_index(read_counts)

        cn_matrix[cell_id] = read_counts['reads']

    cn_matrix = pd.DataFrame(cn_matrix)

    cell_data = pd.DataFrame({'cell_id': cn_matrix.columns.values}).set_index('cell_id')

    adata = ad.AnnData(
        cn_matrix.T,
        obs=cell_data,
        var=bin_data,
    )

    return adata


def read_medicc2_cn(cn_profiles_filename, allele_specific: bool = False) -> AnnData:
    """ Read medicc2 results

    Parameters
    ----------
    cn_profiles_filename : str
        Copy number profiles filename
    allele_specific : bool, optional
        _description_, by default False

    Returns
    -------
    AnnData
        Medicc CN results.
    """

    cn_data = pd.read_csv(
        cn_profiles_filename,
        sep='\t',
        dtype={
            'chrom': 'category',
            'sample_id': 'category'
    })

    cn_data = cn_data.rename(columns={
        'sample_id': 'cell_id',
        'chrom': 'chr',
    })

    if allele_specific:
        cn_data['state'] = cn_data['cn_a'] + cn_data['cn_b']
        cn_fields = ['cn_a', 'cn_b']

    else:
        cn_data = cn_data.rename(columns={'cn': 'state'})
        cn_fields = []

    cn_data['chr'] = cn_data['chr'].str.replace('chr', '')
    cn_data.loc[cn_data['chr'] == '23', 'chr'] = 'X'
    cn_data.loc[cn_data['chr'] == '24', 'chr'] = 'Y'
    cn_data['chr'] = cn_data['chr'].astype('category')

    cn_data['bin'] = cn_data['chr'].astype(str) + ':' + cn_data['start'].astype(str) + '-' + cn_data['end'].astype(str)

    cn_matrix = (
        cn_data
            .set_index(['bin', 'cell_id'])[cn_fields + ['state', 'is_gain', 'is_loss']]
            .unstack(level='cell_id')
            .transpose())

    bin_data = (
        cn_data[['bin', 'chr', 'start', 'end']]
            .drop_duplicates(subset=['bin'])
            .set_index(['bin'])
            .reindex(cn_matrix.loc['state'].columns))

    cell_data = (
        cn_data[['cell_id', 'is_wgd', 'is_normal', 'is_clonal']]
            .drop_duplicates()
            .set_index(['cell_id'])
            .reindex(cn_matrix.loc['state'].index))

    cell_data['is_root'] = cell_data.index == 'diploid'
    cell_data['is_internal'] = cell_data.index.to_series().str.startswith('internal_')
    cell_data['is_cell'] = (~cell_data['is_root']) & (~cell_data['is_internal'])

    layers = {
        'is_gain': cn_matrix.loc['is_gain'],
        'is_loss': cn_matrix.loc['is_loss'],
    }

    for field in cn_fields:
        layers[field] = cn_matrix.loc[field]

    adata = ad.AnnData(
        cn_matrix.loc['state'],
        obs=cell_data,
        var=bin_data,
        layers=layers,
    )

    return adata
