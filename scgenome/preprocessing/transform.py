import numpy as np

from anndata import AnnData


def fill_missing(adata: AnnData, layer_name=None) -> AnnData:
    """ Fill missing values.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default None

    Returns
    -------
    AnnData
        copy number data with no missing data in layer `layer_name`
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

    return adata


