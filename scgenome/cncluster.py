import numpy as np
import scipy.spatial


def compute_bic(kmeans, X):
    """ Computes the BIC metric for a given k means clustering

    Args:
        kmeans: a fitted kmeans clustering object
        X: data for which to calculate bic
    
    Returns:
        float: bic
    
    Reference: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
    """
    centers = [kmeans.cluster_centers_]
    labels = kmeans.labels_
    n_clusters = kmeans.n_clusters
    cluster_sizes = np.bincount(labels)
    N, d = X.shape

    # Compute variance for all clusters
    cl_var = (1.0 / (N - n_clusters) / d) * sum(
        [sum(scipy.spatial.distance.cdist(X[np.where(labels == i)], [centers[0][i]],
                                          'euclidean') ** 2) for i in range(n_clusters)])

    const_term = 0.5 * n_clusters * np.log(N) * (d + 1)

    bic = np.sum([cluster_sizes[i] * np.log(cluster_sizes[i]) -
                  cluster_sizes[i] * np.log(N) -
                  ((cluster_sizes[i] * d) / 2) * np.log(2 * np.pi * cl_var) -
                  ((cluster_sizes[i] - 1) * d / 2) for i in range(n_clusters)]) - const_term

    return bic
