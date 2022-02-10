import sklearn.decomposition
import pandas as pd
import anndata as ad
import numpy as np

import scgenome.preprocessing.transform

from anndata import AnnData


def pca_loadings(adata: AnnData, layer=None, n_components=None, random_state=100) -> AnnData:
    """ Compute PCA loadings matrix

    Parameters
    ----------
    adata : AnnData
        Copy number or other cell matrix
    layer : str, optional
        layer to use, by default None, use .X
    n_components : int, optional
        sklearn.decomposition.PCA n_components parameter, by default None
    random_state : int, optional
        sklearn.decomposition.PCA random_state parameter, by default 100

    Returns
    -------
    DataFrame
        dataframe of pca loadings (rows) by cell (columns)
    """    
    
    """ Calculate loadings of PCA decomposition of a feature matrix.

    Args:
        adata (anndata.AnnData): feature matrix
    """

    if layer is None:
        data = adata.X
    else:
        data = adata.layers[layer]

    pca = sklearn.decomposition.PCA(n_components=n_components, random_state=random_state)
    pca.fit(scgenome.preprocessing.transform.fill_missing(data))

    var = adata.var.copy()
    var['pca_mean'] = pca.mean_

    obs = pd.DataFrame({
        'component': 'PC' + pd.Series(np.arange(pca.explained_variance_.shape[0]) + 1).astype(str),
        'pca_explained_variance': pca.explained_variance_,
        'pca_explained_variance_ratio': pca.explained_variance_ratio_,
        'pca_singular_values': pca.singular_values_,
    }).set_index('component')

    uns = {
        'pca': {
            'params': {
                'layer': layer,
                'n_components': n_components,
                'random_state': random_state,
            },
            'results': {
                'n_components': pca.n_components_,
                'n_features': pca.n_features_,
                'n_samples': pca.n_samples_,
                'noise_variance': pca.noise_variance_,
                'n_features_in': pca.n_features_in_,
            }
        }
    }

    adata = ad.AnnData(
        pca.components_,
        obs=obs,
        var=var,
        uns=uns,
    )

    return adata
