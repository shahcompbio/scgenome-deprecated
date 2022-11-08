import logging
import pandas as pd
import numpy as np
import anndata as ad
import pyranges as pr

import scgenome.loaders.qc
import csverve

from anndata import AnnData
from pyranges import PyRanges
from typing import Dict, Sequence
from pandas import DataFrame


def read_dlp_hmmcopy(alignment_results_dir, hmmcopy_results_dir, annotation_results_dir, sample_ids=None, additional_hmmcopy_reads_cols=None) -> AnnData:
    """ Read hmmcopy results from the DLP pipeline.

    Parameters
    ----------
    alignment_results_dir : str
        dlp pipeline alignment results directory
    hmmcopy_results_dir : str
        dlp pipeline hmmcopy results directory
    annotation_results_dir : str
        dlp pipeline annotation results directory
    sample_ids : list
        sample ids to load
    additional_hmmcopy_reads_cols : list
        per bin metrics to load

    Returns
    -------
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
    ----------
    reads_filename (str):
        dlp pipeline reads filename
    metrics_filename (str):
        dlp pipeline metrics filename
    sample_ids (list):
        sample ids to load

    Returns
    -------
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


def create_cn_anndata(
        cn_data: DataFrame,
        X_column: str,
        layers_columns: Sequence[str],
        cell_metrics_data: DataFrame=None,
        bin_metrics_data: DataFrame=None,
    ) -> AnnData:
    """ Convert hmmcopy pandas dataframes to anndata

    Parameters
    ----------
    cn_data : DataFrame
        copy number data per cell in long format
    X_column : str
        column of cn_data to use for X
    layers_columns : Sequence[str]
        columns of cn_data to use for lauers
    cell_metrics_data : DataFrame, optional
        per cell metrics data, by default None
    bin_metrics_data : DataFrame, optional
        per bin metrics data, by default None

    Returns
    -------
    AnnData
        An instantiated AnnData Object.

    Raises
    ------
    ValueError
        duplicate data or otherwise incompatible inputs
    """
    if cell_metrics_data is None:
        cell_metrics_data = cn_data[['cell_id']].drop_duplicates()

    if bin_metrics_data is None:
        bin_metrics_data = cn_data[['chr', 'start', 'end']].drop_duplicates()

    duplicate_cell_ids = cn_data.loc[cn_data[['chr', 'start', 'end', 'cell_id']].duplicated(keep=False), 'cell_id'].unique()
    if len(duplicate_cell_ids) > 0:
        raise ValueError(f'cell {duplicate_cell_ids[0]} is duplicated, and {len(duplicate_cell_ids)} others')

    assert not cell_metrics_data.duplicated(subset=['cell_id']).any()
    assert not bin_metrics_data.duplicated(subset=['chr', 'start', 'end']).any()

    X = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])[X_column]
            .unstack(level='cell_id')
            .transpose())
    X.index = X.index.astype(str)

    chr_start_end_index = X.columns
    bin_index = (
        X.columns.get_level_values('chr').astype(str) + ':' +
        X.columns.get_level_values('start').astype(str) + '-' +
        X.columns.get_level_values('end').astype(str))

    X = X.set_axis(bin_index, axis=1, copy=False)

    layers = {}
    for layer_name in layers_columns:
        layers[layer_name] = (
            cn_data
                .set_index(['chr', 'start', 'end', 'cell_id'])[layer_name]
                .unstack(level='cell_id')
                .transpose()
                .reindex(index=X.index, columns=chr_start_end_index)
                .set_axis(bin_index, axis=1, copy=False))
        layers[layer_name].index = layers[layer_name].index.astype(str)

    bin_data = (
        bin_metrics_data
            .set_index(['chr', 'start', 'end'], drop=False)
            .reindex(index=chr_start_end_index)
            .set_axis(bin_index, axis=0))

    cell_data = (
        cell_metrics_data
            .set_index(['cell_id'])
            .reindex(X.index))
    cell_data.index = cell_data.index.astype(str)

    adata = ad.AnnData(
        X,
        dtype=X.dtypes[0],
        obs=cell_data,
        var=bin_data,
        layers=layers,
    )

    return adata


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

    return create_cn_anndata(
        cn_data[['cell_id', 'chr', 'start', 'end', 'reads', 'copy', 'state']],
        layers_columns = ['copy', 'state'],
        X_column = 'reads',
        cell_metrics_data=metrics_data,
        bin_metrics_data=cn_data[['chr', 'start', 'end', 'gc']].drop_duplicates())


def convert_dlp_signals(hscn: DataFrame, metrics_data: DataFrame) -> AnnData:
    """ Convert signals pandas dataframes to anndata

    Parameters
    ----------
    hscn : DataFrame
        signals reads data

    Returns
    -------
    AnnData
        An instantiated AnnData Object.
    """
    hscn['state_a'] = hscn['state_AS_phased'].str.split('|', expand=True)[0].astype(float)
    hscn['state_b'] = hscn['state_AS_phased'].str.split('|', expand=True)[1].astype(float)

    layers_columns = [
        'copy', 'state',
        'alleleA', 'alleleB', 'BAF',
        'Maj', 'Min',
        'state_a', 'state_b',
    ]

    return create_cn_anndata(
        hscn,
        layers_columns=layers_columns,
        X_column='totalcounts',
        cell_metrics_data=metrics_data,
    )


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
    bams : Dict[str]
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

        assert (read_counts.index == bin_data.index).all()

        cn_matrix[cell_id] = read_counts['reads'].values.copy()

    cn_matrix = pd.DataFrame(cn_matrix, index=bin_data.index)

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
        cn_data[['bin', 'chr', 'start', 'end', 'is_normal', 'is_clonal']]
            .drop_duplicates(subset=['bin'])
            .set_index(['bin'])
            .reindex(cn_matrix.loc['state'].columns))

    cell_data = (
        cn_data[['cell_id', 'is_wgd']]
            .drop_duplicates()
            .set_index(['cell_id'])
            .reindex(cn_matrix.loc['state'].index))

    cell_data['is_root'] = cell_data.index == 'diploid'
    cell_data['is_internal'] = cell_data.index.to_series().str.startswith('internal_')
    cell_data['is_cell'] = (~cell_data['is_root']) & (~cell_data['is_internal'])

    X = cn_matrix.loc['state'].astype(np.float32)

    layers = {
        'is_gain': cn_matrix.loc['is_gain'],
        'is_loss': cn_matrix.loc['is_loss'],
    }

    for field in cn_fields:
        layers[field] = cn_matrix.loc[field]

    X.index = X.index.astype(str)
    cell_data.index = cell_data.index.astype(str)
    for field in layers:
        layers[field].index = layers[field].index.astype(str)

    adata = ad.AnnData(
        X,
        obs=cell_data,
        var=bin_data,
        layers=layers,
    )

    return adata
