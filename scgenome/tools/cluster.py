import logging
import sklearn.cluster
import sklearn.mixture
import sklearn.preprocessing
import umap
import pandas as pd
import numpy as np
import anndata as ad
from natsort import natsorted

from anndata import AnnData
from typing import Dict, Any, Union
from collections.abc import Iterable

import scgenome.cncluster
import scgenome.preprocessing.transform


def _kmeans_bic(X, k):
    model = sklearn.cluster.KMeans(n_clusters=k, init='k-means++', random_state=100).fit(X)
    bic = scgenome.cncluster.compute_bic(model, X)
    labels = model.labels_

    return labels, bic


def _gmm_diag_bic(X, k):
    model = sklearn.mixture.GaussianMixture(n_components=k, covariance_type='diag', init_params='kmeans')
    labels = model.fit_predict(X)
    bic = model.bic(X)

    return labels, bic


def cluster_cells(
        adata: AnnData,
        layer_name: Union[None, str, Iterable[Union[None,str]]]='copy',
        method: str='kmeans_bic',
        min_k: int=2,
        max_k: int=100,
        cell_ids: Iterable[str]=None,
        bin_ids: Iterable[str]=None,
        standardize: bool=False,
    ) -> AnnData:
    """ Cluster cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    method : str, optional
        clustering method, by default 'kmeans_bic'
    min_k : int, optional
        minimum number of clusters, by default 2
    max_k : int, optional
        maximum number of clusters, by default 100
    cell_ids : str, optional
        subset of cells to cluster, by default None
    bin_ids : str, optional
        subset of bins to cluster, by default None
    standarize : bool
        standardize the data prior to outlier detection, by default False
    Returns
    -------
    AnnData
        copy number data with additional `cluster_id` and `cluster_size` columns

    Examples
    -------

    >>> import scgenome
    >>> import anndata as ad
    >>> import numpy as np
    >>> adata = ad.AnnData(np.array([
    ...    [3, 3, 3, 6, 6],
    ...    [1, 1, 1, 2, 2],
    ...    [1, 22, 1, 2, 2],
    ...    [1, 3, 3, 5, 5],
    ... ]).astype(np.float32))
    >>> adata = scgenome.tl.cluster_cells_kmeans(adata, layer_name=None, max_k=3)
    >>> adata.obs['cluster_id']
    0    0
    1    2
    2    1
    3    0
    Name: cluster_id, dtype: category
    Categories (3, int64): [0, 1, 2]

    """
    if cell_ids is None:
        cell_ids = adata.obs.index

    if bin_ids is None:
        bin_ids = adata.var.index

    min_k = min(adata.shape[0], min_k)
    max_k = min(adata.shape[0], max_k)

    def __get_layer(layer_name):
        if layer_name is not None:
            return np.array(adata[cell_ids, bin_ids].layers[layer_name])
        else:
            return np.array(adata[cell_ids, bin_ids].X)

    if isinstance(layer_name, str):
        X = __get_layer(layer_name)
    elif isinstance(layer_name, Iterable):
        X = np.concatenate([__get_layer(l) for l in layer_name], axis=1)

    X = scgenome.preprocessing.transform.fill_missing(X)
    
    if standardize:
        X = sklearn.preprocessing.StandardScaler().fit_transform(X)

    ks = range(min_k, max_k + 1)

    logging.info(f'trying with max k={max_k}')

    labels = []
    criterias = []
    for k in ks:
        logging.info(f'trying with k={k}')
        if method == 'kmeans_bic':
            label, criteria = _kmeans_bic(X, k)
        elif method == 'gmm_diag_bic':
            label, criteria = _gmm_diag_bic(X, k)
        else:
            raise ValueError(f'unrecognized method {method}')
        labels.append(label)
        criterias.append(criteria)

    opt_k_idx = np.array(criterias).argmax()
    opt_k = ks[opt_k_idx]
    opt_label = labels[opt_k_idx]
    logging.info(f'selected k={opt_k}')
    
    adata.obs['cluster_id'] = '-1'
    adata.obs.loc[cell_ids, 'cluster_id'] = pd.Series(opt_label, index=adata.obs.loc[cell_ids].index).astype('str').astype('category')
    adata.obs['cluster_size'] = adata.obs.groupby('cluster_id')['cluster_id'].transform('size')

    # store information on the clustering parameters
    adata.uns['clustering'] = {}
    adata.uns['clustering']['params'] = dict(
        method=method,
        opt_k=opt_k,
        min_k=min_k,
        max_k=max_k,
        layer_name=layer_name,
        cell_ids=np.array(cell_ids),
        bin_ids=np.array(bin_ids),
        standardize=standardize,
    )

    return adata


def detect_outliers(
        adata: AnnData,
        layer_name: Union[None, str, Iterable[Union[None,str]]]='copy',
        method: str='isolation_forest',
        standarize: bool=False,
    ) -> AnnData:
    """ Detect outlier cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to use for outlier detection, None for X, by default 'copy'
    method : str, optional
        outlier method, by default 'isolation_forest'
    standarize : bool
        standardize the data prior to outlier detection, by default False

    Returns
    -------
    AnnData
        copy number data with additional `is_outlier` column

    """
    def __get_layer(layer_name):
        if layer_name is not None:
            return np.array(adata.layers[layer_name])
        else:
            return np.array(adata.X)

    if isinstance(layer_name, str):
        X = __get_layer(layer_name)
    elif isinstance(layer_name, Iterable):
        X = np.concatenate([__get_layer(l) for l in layer_name], axis=1)

    X = scgenome.preprocessing.transform.fill_missing(X)

    if standarize:
        X = sklearn.preprocessing.StandardScaler().fit_transform(X)

    if method == 'isolation_forest':
        model = sklearn.ensemble.IsolationForest()
        is_outlier = (model.fit_predict(X) == -1) * 1
    elif method == 'local_outlier_factor':
        model = sklearn.neighbors.LocalOutlierFactor()
        is_outlier = (model.fit_predict(X) == -1) * 1
    else:
        raise ValueError(f'unknown method {method}')

    adata.obs['is_outlier'] = pd.Series(is_outlier, index=adata.obs.index, dtype='category')

    # store information on the clustering parameters
    adata.uns['outliers'] = {}
    adata.uns['outliers']['params'] = dict(
        method='isolation_forest',
        layer_name=layer_name,
        standarize=standarize,
    )

    return adata


def aggregate_clusters(
        adata: AnnData,
        agg_X: Any,
        agg_layers: Dict=None,
        agg_obs: Dict=None,
        cluster_col: str='cluster_id',
        cluster_size_col: str='cluster_size') -> AnnData:
    """ Aggregate copy number by cluster to create cluster CN matrix

    Parameters
    ----------
    adata : AnnData
        copy number data
    agg_X : Any
        function to aggregate X
    agg_layers : Dict, optional
        functions to aggregate layers keyed by layer names, by default None
    agg_obs : Dict, optional
        functions to aggregate obs data keyed by obs columns, by default None
    cluster_col : str, optional
        column with cluster ids, by default 'cluster_id'
    cluster_size_col : str, optional
        column that will be set to the size of each cluster, by default 'cluster_size'

    Returns
    -------
    AnnData
        aggregated cluster copy number
    """

    X = (
        adata
            .to_df()
            .set_index(adata.obs[cluster_col].astype(str))
            .groupby(level=0)
            .agg(agg_X)
            .sort_index())

    layer_data = None
    if agg_layers is not None:
        layer_data = {}
        for layer_name in agg_layers:
            layer_data[layer_name] = (
                adata
                    .to_df(layer=layer_name)
                    .set_index(adata.obs[cluster_col].astype(str))
                    .groupby(level=0)
                    .agg(agg_layers[layer_name])
                    .sort_index())

    obs_data = {}
    obs_data[cluster_size_col] = (
        adata.obs
            .set_index(adata.obs[cluster_col].astype(str))
            .groupby(level=0)
            .size())

    if agg_obs is not None:
        for obs_name in agg_obs:
            obs_data[obs_name] = (
                adata.obs
                    .set_index(adata.obs[cluster_col].astype(str))[obs_name]
                    .groupby(level=0)
                    .agg(agg_obs[obs_name])
                    .sort_index())

    obs_data = pd.DataFrame(obs_data)

    dtypes = X.dtypes.unique()
    assert len(dtypes) == 1
    dtype = dtypes[0]

    adata = ad.AnnData(
        X,
        dtype=dtype,
        obs=obs_data,
        var=adata.var,
        layers=layer_data,
    )

    return adata


def aggregate_clusters_hmmcopy(adata: AnnData) -> AnnData:
    """ Aggregate hmmcopy copy number by cluster to create cluster CN matrix

    Parameters
    ----------
    adata : AnnData
        hmmcopy copy number data

    Returns
    -------
    AnnData
        aggregsated cluster copy number
    """

    agg_X = np.sum

    agg_layers = {
        'copy': np.nanmean,
        'state': np.nanmedian,
    }

    agg_obs = {
        'total_reads': np.nansum,
    }

    return aggregate_clusters(adata, agg_X, agg_layers, agg_obs, cluster_col='cluster_id')


def compute_umap(
        adata: AnnData,
        layer_name: str='copy',
        n_components: int=2,
        n_neighbors: int=15,
        min_dist: float=0.1,
        metric: str='euclidean',
    ) -> AnnData:
    """ Cluster cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data on which to perform umap dimensionality
        reduction, None for X, by default 'copy'
    n_components : int
        umap n_components param
    n_neighbors : int
        umap n_neighbors param
    min_dist : float
        umap min_dist param

    Returns
    -------
    AnnData
        copy number data with additional `umap_1`, `umap_2` columns
    """

    if layer_name is not None:
        X = adata.layers[layer_name]
    else:
        X = adata.X

    X = scgenome.preprocessing.transform.fill_missing(X)

    embedding = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=42,
    ).fit_transform(X)

    adata.obs['UMAP1'] = embedding[:, 0]
    adata.obs['UMAP2'] = embedding[:, 1]

    return adata
