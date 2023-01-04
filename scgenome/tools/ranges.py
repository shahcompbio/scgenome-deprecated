import pyranges as pr
import pandas as pd
import anndata as ad

from pyranges import PyRanges
from pandas import DataFrame
from anndata import AnnData
from numpy import ndarray

import scgenome.refgenome



def dataframe_to_pyranges(data: DataFrame) -> PyRanges:
    """ Convert dataframe with chr, start, end to pyranges object

    Parameters
    ----------
    data : DataFrame
        dataframe with columns 'chr', 'start', 'end'

    Returns
    -------
    PyRanges
        pyranges object
    """
    data = pr.PyRanges(data.rename(columns={
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End',
    }))

    return data


def pyranges_to_dataframe(data: PyRanges) -> DataFrame:
    """ Convert pyranges object to dataframe

    Parameters
    ----------
    data : PyRanges
        pyranges object

    Returns
    -------
    DataFrame
        dataframe with columns 'chr', 'start', 'end'
    """
    data = data.as_df().rename(columns={
        'Chromosome': 'chr',
        'Start': 'start',
        'End': 'end',
    })

    return data


def create_bins_pr(binsize: int) -> PyRanges:
    """ Create a regular binning of the genome

    Parameters
    ----------
    binsize : int
        size of bins

    Returns
    -------
    PyRanges
        regular bins tiled across the genome
    """    
    
    chromsizes = scgenome.refgenome.info.chromosome_info[['chr', 'chromosome_length']].rename(
        columns={'chr': 'Chromosome', 'chromosome_length': 'End'}).assign(Start=0)[['Chromosome', 'Start', 'End']]

    chromsizes = pr.PyRanges(chromsizes)

    bins = pr.gf.tile_genome(chromsizes, binsize)

    return bins


def create_bins(binsize: int) -> DataFrame:
    """ Create a regular binning of the genome

    Parameters
    ----------
    binsize : int
        size of bins

    Returns
    -------
    DataFrame
        regular bins tiled across the genome
    """
    return pyranges_to_dataframe(create_bins_pr(binsize=binsize))


def rebin_agg_df(data: DataFrame, intersect: DataFrame, agg_f: dict) -> DataFrame:
    """ Rebin a dataframe, aggregating across multiple intersecting regions

    Parameters
    ----------
    data : DataFrame
        binned data to rebin and aggregate with columns 'chr', 'start', 'end' and index 'bin'
    intersect : DataFrame
        mapping between previous and target bins with columns 'bin', 'target_bin', 'width'
    agg_f : dict of tuples
        aggregate functions, similar to pandas dataframe groupby .agg

    Returns
    -------
    DataFrame
        rebinned and aggregated data with columns according to 'agg_f' and index 'target_bin'
    """
    data = data.loc[intersect['bin'], :]

    data = data.set_index(intersect.set_index('bin')['target_bin'], append=True)
    data = data.set_index(intersect.set_index('bin')['width'], append=True)
    
    data = data.groupby(level='target_bin').agg(**agg_f)

    return data


def rebin_agg_layer(adata: AnnData, intersect: DataFrame, layer_name: str, agg_f: dict) -> DataFrame:
    """ Rebin an anndata layer, aggregating across multiple intersecting regions

    Parameters
    ----------
    adata : DataFrame
        adata to rebin and aggregate with var columns 'chr', 'start', 'end' and index 'bin'
    intersect : DataFrame
        mapping between previous and target bins with columns 'bin', 'target_bin', 'width'
    layer_name : str
        name of layer to rebin and aggregate
    agg_f : dict of tuples
        aggregate functions, similar to pandas dataframe groupby .agg

    Returns
    -------
    DataFrame
        rebinned and aggregated layer data with columns according to 'agg_f' and index 'target_bin'
    """
    data = adata.to_df(layer=layer_name)

    data = data.loc[:, intersect['bin']].T

    data = data.set_index(intersect.set_index('bin')['target_bin'], append=True)
    data = data.set_index(intersect.set_index('bin')['width'], append=True)

    data = data.groupby(level='target_bin').agg(agg_f)
    
    return data.T


def intersect_regions(a: DataFrame, b: DataFrame) -> DataFrame:
    """ Compute intersection between two sets of regions

    Parameters
    ----------
    a : DataFrame
        region data with columns 'chr', 'start', 'end', ...
    b : DataFrame
        region data with columns 'chr', 'start', 'end', ...

    Returns
    -------
    DataFrame
        intersecting regions with columns 'chr', 'start', 'end' and columns from 'a' and 'b'
    """
    a = scgenome.tools.ranges.dataframe_to_pyranges(a)
    b = scgenome.tools.ranges.dataframe_to_pyranges(b)

    intersect_1 = a.intersect(b)
    intersect_2 = b.intersect(a)

    intersect = pd.merge(
        scgenome.tools.ranges.pyranges_to_dataframe(intersect_1),
        scgenome.tools.ranges.pyranges_to_dataframe(intersect_2),
        on=['chr', 'start', 'end'])

    return intersect


def rebin(adata: AnnData, target_bins: DataFrame, outer_join: bool=False, agg_X=None, agg_var=None, agg_layers=()) -> AnnData:
    """ Rebin an AnnData and aggregate across multiple intersecting regions

    Parameters
    ----------
    adata : AnnData
        data to rebin and aggregate
    target_bins : DataFrame
        target bins with columns 'chr', 'start', 'end'
    outer_join : bool, optional
        whether to include target bins with no overlap in adata, by default False
    agg_X : callable, optional
        aggregate function for X, by default None
    agg_var : dict, optional
        aggregate functions for var, by default None
    agg_layers : dict, optional
        aggregate functions for each layer, by default ()

    Returns
    -------
    AnnData
        rebinned and aggregated data
    """
    bins = adata.var.rename_axis('bin').reset_index()[['chr', 'start', 'end', 'bin']]

    target_bins['target_bin'] = (
        target_bins['chr'].astype(str) + ':' +
        target_bins['start'].astype(str) + '-' +
        target_bins['end'].astype(str))
    target_bins = target_bins[['chr', 'start', 'end', 'target_bin']]

    intersect = intersect_regions(bins, target_bins)
    intersect['width'] = intersect['end'] - intersect['start'] + 1

    if agg_var is not None:
        var = rebin_agg_df(adata.var, intersect, agg_var)
    
    else:
        var = intersect[['target_bin']].set_index('target_bin')

    var = target_bins.set_index('target_bin').merge(
        var, left_index=True, right_index=True, how=('inner', 'left')[outer_join])

    var = var.sort_values(['chr', 'start'])

    if agg_X is not None:
        X = rebin_agg_layer(adata, intersect, None, agg_X)
        X = X.reindex(columns=var.index)

    else:
        X = None

    layer_data = {}
    for layer_name in agg_layers:
        layer_data[layer_name] = rebin_agg_layer(adata, intersect, layer_name, agg_layers[layer_name])
        layer_data[layer_name] = layer_data[layer_name].reindex(columns=var.index)

    adata = ad.AnnData(
        X,
        obs=adata.obs,
        var=var,
        layers=layer_data,
    )

    return adata


def rebin_regular(adata: AnnData, bin_size: int, outer_join: bool=False, agg_X=None, agg_var=None, agg_layers=()) -> AnnData:
    """ Rebin an AnnData and aggregate across multiple intersecting regions

    Parameters
    ----------
    adata : AnnData
        data to rebin and aggregate
    bin_size : int
        width of target bins
    outer_join : bool, optional
        whether to include target bins with no overlap in adata, by default False
    agg_X : callable, optional
        aggregate function for X, by default None
    agg_var : dict, optional
        aggregate functions for var, by default None
    agg_layers : dict, optional
        aggregate functions for each layer, by default ()

    Returns
    -------
    AnnData
        rebinned and aggregated data
    """
    target_bins = create_bins(bin_size)
    target_bins['start'] = target_bins['start'] + 1

    return rebin(adata, target_bins, outer_join=outer_join, agg_X=agg_X, agg_var=agg_var, agg_layers=agg_layers)


def weighted_mean(values: ndarray, widths: ndarray) -> float:
    """ Compute weighted mean

    Parameters
    ----------
    values : ndarray
        values to compute mean 
    widths : ndarray
        weights corresponding to each value

    Returns
    -------
    float
        weighted mean
    """
    weighted_mean = (widths * values).sum() / widths.sum()
    return weighted_mean


def bin_width_weighted_mean(df: DataFrame) -> float:
    """ Convenience function to allow computing of weighted mean with rebin aggregate functions

    Parameters
    ----------
    df : DataFrame
        group of data for which to compute mean, must have width in index

    Returns
    -------
    float
        weighted mean
    """
    values = df.values
    widths = df.index.get_level_values('width').values
    return weighted_mean(values, widths)


