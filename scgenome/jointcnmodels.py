import numpy as np
from remixt.bpmodel import sum_product
from remixt.bpmodel import sum_product_2paramtrans
import scipy


def calculate_ll_normal_simple(data, variances):
    """ Calculate likelihood per state per segment
    
    Args:
        data: (n_cells, n_segments)
        variances: (n_cells, n_states)

    Returns:
        log likelihood: (n_cells, n_segments, n_states)
    """

    n_cells = data.shape[0]
    n_segments = data.shape[1]
    n_states = variances.shape[1]

    # Create states as (n_segments, n_cells, n_states) array
    states = np.tile(np.arange(0, n_states, 1), (n_cells, n_segments, 1))

    # Calculate mean
    mean = states
    
    # Normal dist log likelihood
    ll = (
        -0.5 * np.log(2. * np.pi)
        -0.5 * np.log(variances[:, np.newaxis, :])
        -1. * (np.square(data[:, :, np.newaxis] - mean) /
                (2. * variances[:, np.newaxis, :])))
    ll[np.isnan(data)] = 0.

    return ll


def calculate_marginal_ll_simple(data, variances, transmodel):
    framelogprob = calculate_ll_normal_simple(data, variances).sum(axis=0)

    alphas = np.zeros(framelogprob.shape)
    betas = np.zeros(framelogprob.shape)

    if transmodel['kind'] == 'twoparam':
        sum_product_2paramtrans(
            framelogprob,
            alphas,
            betas,
            transmodel['e0'],
            transmodel['e1'],
        )

    elif transmodel['kind'] == 'full':
        sum_product(
            framelogprob,
            transmodel['tr_mat'],
            alphas,
            betas,
        )

    else:
        raise ValueError("unknown transition model {transmodel['kind']}")

    return scipy.special.logsumexp(alphas[-1, :])


def gibbs_sample_cluster_indices(data, variances, assignments, max_clusters, alpha, transmodel):
    n_cells = data.shape[0]

    for cell_idx in range(n_cells):
        log_assign_prob = np.zeros((max_clusters,))

        for cluster_idx in range(max_clusters):
            assignments[cell_idx] = cluster_idx
            log_marginal_with = calculate_marginal_ll_simple(
                data[assignments == cluster_idx, :],
                  variances[assignments == cluster_idx, :],
                transmodel)

            assignments[cell_idx] = -1
            log_marginal_without = calculate_marginal_ll_simple(
                data[assignments == cluster_idx, :],
                variances[assignments == cluster_idx, :],
                transmodel)

            log_posterior_predictive = log_marginal_with - log_marginal_without
            num_cluster = np.sum(assignments[cell_idx] == cluster_idx)
            if num_cluster == 0:
                log_assign_prob[cluster_idx] = np.log(alpha / (alpha + n_cells - 1)) + log_posterior_predictive
            else:
                log_assign_prob[cluster_idx] = np.log(num_cluster / (alpha + n_cells - 1)) + log_posterior_predictive
        
        assign_prob = np.exp(log_assign_prob - log_assign_prob.max())
        assign_prob /= assign_prob.sum()

        assignments[cell_idx] = np.random.choice(assign_prob.shape[0], p=assign_prob)
    
    return assignments


