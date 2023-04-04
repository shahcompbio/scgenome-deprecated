import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dst
import sklearn.preprocessing
import pandas as pd

from anndata import AnnData
from collections.abc import Iterable
from typing import Union, Any, Dict

import scgenome.preprocessing.transform


def sort_cells(
        adata: AnnData,
        layer_name: Union[None, str, Iterable[Union[None,str]]]='copy',
        cell_ids: Iterable[str]=None,
        bin_ids: Iterable[str]=None,
        standarize: bool=False,
    ) -> AnnData:
    """ Sort cells by hierarchical clustering on copy number values.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to use for sorting, None for X, by default 'copy'
    cell_ids : str, optional
        subset of cells to cluster, by default None
    bin_ids : str, optional
        subset of bins to cluster, by default None
    standarize : bool
        standardize the data prior to sorting, by default False

    Returns
    -------
    AnnData
        copy number data with cell_order column added to obs
    """
    if cell_ids is None:
        cell_ids = adata.obs.index

    if bin_ids is None:
        bin_ids = adata.var.index

    def __get_layer(layer_name):
        if layer_name is not None:
            return np.array(adata[cell_ids, bin_ids].layers[layer_name])
        else:
            return np.array(adata[cell_ids, bin_ids].X)

    if isinstance(layer_name, (str, type(None))):
        X = __get_layer(layer_name)
    elif isinstance(layer_name, Iterable):
        X = np.concatenate([__get_layer(l) for l in layer_name], axis=1)
    else:
        raise ValueError(f'layer_name was {layer_name}')

    X = scgenome.preprocessing.transform.fill_missing(X)

    if standarize:
        X = sklearn.preprocessing.StandardScaler().fit_transform(X)

    D = dst.squareform(dst.pdist(X, 'cityblock'))
    Y = sch.linkage(D, method='complete')
    Z = sch.dendrogram(Y, color_threshold=-1, no_plot=True)
    idx = np.array(Z['leaves'])

    ordering = np.zeros(idx.shape[0], dtype=int)
    ordering[idx] = np.arange(idx.shape[0])

    adata.obs['cell_order'] = '-1'
    adata.obs.loc[cell_ids, 'cell_order'] = pd.Series(ordering, index=adata.obs.loc[cell_ids].index)

    return adata


def sort_clusters(
        adata: AnnData,
        layer_name: Union[None, str, Iterable[Union[None,str]]]='copy',
        agg_X: Any=None,
        agg_layers: Dict=None,
        cell_ids: Iterable[str]=None,
        bin_ids: Iterable[str]=None,
        standarize: bool=False,
    ) -> AnnData:
    """ Sort cells by hierarchical clustering on copy number values.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to use for sorting, None for X, by default 'copy'
    agg_X : Any
        function to aggregate X, by default None
    agg_layers : Dict, optional
        functions to aggregate layers keyed by layer names, by default None
    cell_ids : str, optional
        subset of cells to cluster, by default None
    bin_ids : str, optional
        subset of bins to cluster, by default None
    standarize : bool
        standardize the data prior to sorting, by default False

    Returns
    -------
    AnnData
        copy number data with cell_order column added to obs
    """

    if cell_ids is None:
        cell_ids = adata.obs.index

    if bin_ids is None:
        bin_ids = adata.var.index

    adata_clusters = scgenome.tools.cluster.aggregate_clusters(
        adata[cell_ids, bin_ids], agg_X=agg_X, agg_layers=agg_layers)

    adata_clusters = sort_cells(
        adata_clusters,
        layer_name=layer_name,
        standarize=standarize)

    adata.obs['cell_order'] = '-1'
    adata.obs.loc[cell_ids, 'cluster_order'] = pd.Series(
        adata_clusters.obs.loc[adata.obs.loc[cell_ids, 'cluster_id'].values, 'cell_order'].values,
        index=adata.obs.loc[cell_ids].index)

    return adata
