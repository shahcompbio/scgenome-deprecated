import logging
import sklearn.cluster
import umap
import pandas as pd
import numpy as np
from natsort import natsorted

from anndata import AnnData
from pandas import DataFrame

import scgenome.cncluster


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
    """

    if layer_name is not None:
        X = adata.layers[layer_name]
    else:
        X = adata.X

    # Deal with missing values by assigning the mean value
    # of each bin to missing values of that bin
    bin_means = np.nanmean(X, axis=0)
    bin_means = np.nan_to_num(bin_means, nan=0)
    bin_means = np.tile(bin_means, (X.shape[0], 1))
    X[np.where(np.isnan(X))] = bin_means[np.where(np.isnan(X))]

    ks = range(min_k, max_k + 1)

    print(f'trying with max k={max_k}')

    kmeans = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.cluster.KMeans(n_clusters=k, init="k-means++").fit(X)
        bic = scgenome.cncluster.compute_bic(model, X)
        kmeans.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = kmeans[opt_k]

    adata.obs['cluster_id'] = model.labels_

    # store information on the clustering parameters
    adata.uns['kmeans'] = {}
    adata.uns['kmeans']['params'] = dict(
        opt_k=opt_k,
        min_k=min_k,
        max_k=max_k,
        layer_name=layer_name,
    )

    return adata
