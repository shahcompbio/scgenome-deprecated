import pyranges as pr
import pandas as pd
import anndata as ad
import numpy as np

from pandas import DataFrame
from anndata import AnnData
from pyranges import PyRanges
from collections.abc import Iterable


def read_ensemble_genes_gtf(gtf_filename) -> PyRanges:
    """ Read an ensembl gtf and extract gene start end

    Parameters
    ----------
    gtf_filename : str
        GTF filename

    Returns
    -------
    PyRanges
        Genes bounds
    """    
    genes = pr.read_gtf(gtf_filename, as_df=True)
    genes = genes.groupby(['Chromosome', 'gene_id', 'gene_name'], observed=True).agg({'Start': min, 'End': max}).reset_index()
    genes = pr.PyRanges(genes)
    return genes


def _segment_width_weighted_mean_matrix(data: DataFrame, intersect: DataFrame) -> DataFrame:
    # Add gene id to index
    data = data.loc[:, intersect['bin']].T.set_index(pd.Index(intersect['gene_id'].values, name='gene_id'), append=True)

    # Scale by overlap segment width
    data = data * intersect['segment_width'].values[:, np.newaxis]

    # Calculate segment width normalized mean
    data = data.groupby(level=1).sum()
    data = data.T / intersect.groupby('gene_id')['segment_width'].sum()

    return data


def _segment_width_weighted_mean_var(data: DataFrame, intersect: DataFrame) -> DataFrame:
    # Add gene id to index
    data = data.loc[intersect['bin']].set_index(pd.Index(intersect['gene_id'].values, name='gene_id'), append=True)

    # Scale by overlap segment width
    data = data * intersect['segment_width'].values[:, np.newaxis]

    # Calculate segment width normalized mean
    data = data.groupby(level=1).sum()
    data = (data.T / intersect.groupby('gene_id')['segment_width'].sum()).T

    return data


def aggregate_genes(
        adata: AnnData,
        genes: PyRanges,
        agg_layers: Iterable=None,
        agg_var: Iterable=None) -> AnnData:
    """ Aggregate copy number by gene to create gene CN matrix

    Currently only does segment width weighted mean aggregation.

    Parameters
    ----------
    adata : AnnData
        copy number data
    genes : PyRanges
        gene data
    agg_layers : List, optional
        list of layers to aggregate, by default None, all layers
    agg_var : List, optional
        list of obs columns to aggregate, by default None, all columns
    cluster_col : str, optional
        column with cluster ids, by default 'cluster_id'

    Returns
    -------
    AnnData
        aggregated gene copy number
    """

    if agg_layers is None:
        agg_layers = adata.layers.keys()
    agg_layers = set(agg_layers)

    if agg_var is None:
        agg_var = set(adata.var.select_dtypes(include=np.number).columns.to_list()) - set(['chr', 'start', 'end'])
    agg_var = set(agg_var)

    bins = pr.PyRanges(adata.var.reset_index().rename(columns={
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End',
    })[['Chromosome', 'Start', 'End', 'bin']])

    intersect_1 = genes.intersect(bins)
    intersect_2 = bins.intersect(genes)

    intersect = pd.merge(intersect_1.as_df(), intersect_2.as_df())
    intersect['segment_width'] = intersect['End'] - intersect['Start']

    X = _segment_width_weighted_mean_matrix(adata.to_df(), intersect)

    layer_data = {}
    for layer_name in agg_layers:
        layer_data[layer_name] = _segment_width_weighted_mean_matrix(adata.to_df(layer=layer_name), intersect)

    var = _segment_width_weighted_mean_var(adata.var[agg_var], intersect)

    gene_data = genes.as_df().drop_duplicates().set_index('gene_id')
    var = var.merge(gene_data, left_index=True, right_index=True, how='left')

    adata = ad.AnnData(
        X,
        obs=adata.obs,
        var=var,
        layers=layer_data,
    )

    return adata


