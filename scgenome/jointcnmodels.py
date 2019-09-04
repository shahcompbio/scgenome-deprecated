import numpy as np
from remixt.bpmodel import sum_product
import scipy

from scgenome.constants import MAX_CN


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


def calculate_marginal_ll_simple(data, variances, tr_mat):
    # n_bin x n_states
    framelogprob = calculate_ll_normal_simple(data, variances).sum(axis=0)

    alphas = np.zeros(framelogprob.shape)
    betas = np.zeros(framelogprob.shape)

    sum_product(
        framelogprob,
        tr_mat,
        alphas,
        betas)

    return scipy.special.logsumexp(alphas[-1, :])


def gibbs_sample_cluster_indices(data, variances, tr_mat, assignments, max_clusters, alpha):
    n_cells = data.shape[0]

    for cell_idx in range(n_cells):
        log_assign_prob = np.zeros((max_clusters,))

        for cluster_idx in range(max_clusters):
            assignments[cell_idx] = cluster_idx
            log_marginal_with = calculate_marginal_ll_simple(
                data[assignments == cluster_idx, :],
                variances[assignments == cluster_idx, :],
                tr_mat)

            assignments[cell_idx] = -1
            log_marginal_without = calculate_marginal_ll_simple(
                data[assignments == cluster_idx, :],
                variances[assignments == cluster_idx, :],
                tr_mat)

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


def get_variances(cn_data, matrix_data, n_states=MAX_CN):
    cell_state_var = cn_data[['cell_id', 'state', 'copy']].dropna().groupby(
        ['cell_id', 'state'])['copy'].var().rename('copy_var').reset_index()
    variances = cell_state_var.set_index(['state', 'cell_id'])[
        'copy_var'].unstack()
    #variances = variances.reindex(columns=matrix_data['reads'].columns,
    print(matrix_data.shape)
    variances = variances.reindex(columns=matrix_data['copy'].columns,
                                  index=range(n_states)).fillna(0.05).T
    variances = variances.values
    variances[variances < 0.001] = 0.001
    return variances


def get_tr_probs(n_segments, n_states=MAX_CN):
    tr_probs = np.zeros((n_segments, n_states, n_states))
    tr_probs[:] += 1.
    tr_probs[:, range(n_states), range(n_states)] += 100.
    tr_probs /= tr_probs.sum(axis=2)[:, :, np.newaxis]
    return tr_probs
