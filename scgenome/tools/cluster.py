import logging
import sklearn.cluster
import umap
import pandas as pd
import numpy as np
import anndata as ad
from natsort import natsorted

from anndata import AnnData
from typing import Dict, Any

import scgenome.cncluster
import scgenome.preprocessing.transform


def cluster_cells(adata: AnnData, layer_name='copy', method='kmeans') -> AnnData:
    """ Cluster cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    method : str, optional
        clustering method, by default 'kmeans'

    Returns
    -------
    AnnData
        copy number data with additional `cluster_id` column
    """

    if method == 'kmeans':
        cluster_cells_kmeans(adata, layer_name=layer_name)


def cluster_cells_kmeans(adata: AnnData, layer_name='copy', min_k=2, max_k=100) -> AnnData:
    """ Cluster cells by copy number using kmeans.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    min_k : int, optional
        minimum number of clusters, by default 2
    max_k : int, optional
        maximum number of clusters, by default 100

    Returns
    -------
    AnnData
        copy number data with additional `cluster_id` column

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

    if layer_name is not None:
        X = adata.layers[layer_name]
    else:
        X = adata.X

    X = scgenome.preprocessing.transform.fill_missing(X)

    ks = range(min_k, max_k + 1)

    logging.info(f'trying with max k={max_k}')

    kmeans = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.cluster.KMeans(n_clusters=k, init='k-means++', random_state=100).fit(X)
        bic = scgenome.cncluster.compute_bic(model, X)
        kmeans.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = kmeans[opt_k]

    adata.obs['cluster_id'] = pd.Series(model.labels_, index=adata.obs.index, dtype='category')

    # store information on the clustering parameters
    adata.uns['kmeans'] = {}
    adata.uns['kmeans']['params'] = dict(
        opt_k=opt_k,
        min_k=min_k,
        max_k=max_k,
        layer_name=layer_name,
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

    adata = ad.AnnData(
        X,
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
