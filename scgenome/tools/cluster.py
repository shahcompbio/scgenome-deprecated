

def cluster_cells(adata, method='kmeans'):
    """ Cluster cells by copy number.

    Args:
        adata (anndata.AnnData): copy number

    KwArgs:
        method (str): clustering method

    """
    # cncluster.umap_hdbscan_cluster, cncluster.umap_hdbscan_cluster, cnclones.calculate_*_clones

