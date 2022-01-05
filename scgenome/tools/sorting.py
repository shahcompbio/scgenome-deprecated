import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dst

from anndata import AnnData

import scgenome.preprocessing.transform


def sort_cells(adata: AnnData, layer_name=None) -> AnnData:
    """ Sort cells by hierarchical clustering on copy number values.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default None

    Returns
    -------
    AnnData
        copy number data with cell_order column added to obs
    """

    if layer_name is not None:
        X = adata.layers[layer_name]
    else:
        X = adata.X

    X = scgenome.preprocessing.transform.fill_missing(X)

    D = dst.squareform(dst.pdist(X, 'cityblock'))
    Y = sch.linkage(D, method='complete')
    Z = sch.dendrogram(Y, color_threshold=-1, no_plot=True)
    idx = np.array(Z['leaves'])

    ordering = np.zeros(idx.shape[0], dtype=int)
    ordering[idx] = np.arange(idx.shape[0])

    adata.obs['cell_order'] = ordering

    return adata

