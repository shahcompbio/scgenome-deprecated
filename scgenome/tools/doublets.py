import anndata as ad
import pandas as pd
import numpy as np


from .cluster import compute_umap, cluster_cells_kmeans

def infer_doublets(adata, layer_name='copy', synthetic_multiplier=1., fraction_synthetic_threshold=0.95):
    n_cells = adata.shape[0]
    n_synthetic = int(n_cells * synthetic_multiplier)

    pairs = np.array([
        np.random.choice(adata.shape[0], size=2, replace=False) for _ in range(n_synthetic)])

    norm_copy = (adata.layers[layer_name] / 
        np.nansum(adata.layers[layer_name], axis=1)[:, np.newaxis])

    synthetic_doublets = norm_copy[pairs, ].mean(axis=1)

    synthetic_doublets_obs = pd.DataFrame(index=(f'synthetic_{a}' for a in range(synthetic_doublets.shape[0])))

    doublet_adata = ad.AnnData(
        np.concatenate([
            norm_copy,
            synthetic_doublets,
        ]),
        obs=pd.concat([
            adata.obs.assign(synthetic=0.),
            synthetic_doublets_obs.assign(synthetic=1.),
        ]),
        var=adata.var,
    )

    doublet_adata = compute_umap(
        doublet_adata, layer_name=None, metric='correlation')

    doublet_adata = cluster_cells_kmeans(doublet_adata, max_k=30, layer_name=None)

    doublet_adata.obs['fraction_synthetic'] = doublet_adata.obs.groupby('cluster_id')['synthetic'].transform('mean')

    doublet_adata.obs['is_doublet'] = (
        doublet_adata.obs['fraction_synthetic'] > fraction_synthetic_threshold)

    return doublet_adata

