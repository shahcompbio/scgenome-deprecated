import numpy as np

def get_gaussian_bicluster(samples_per_cluster, num_bin):
    mean1 = np.repeat(-1, samples_per_cluster)
    mean2 = np.repeat(1, samples_per_cluster)
    cov = np.diag(samples_per_cluster)
    np.fill_diagonal(cov, 1)
    cluster1 = np.random.multivariate_normal(mean1, cov, num_bin)
    cluster2 = np.random.multivariate_normal(mean2, cov, num_bin)
    cn_mat = np.concatenate([cluster1, cluster2], axis=1).T
    return cn_mat


