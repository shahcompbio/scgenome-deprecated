from scipy.misc import logsumexp as log_sum_exp

import numpy as np
import scipy.stats as stats


def compute_log_likelihoods(df, error_rate=1e-3):
    """ Compute the presence absence log likelihood of an SNV
    """
    df['log_likelihood_absent'] = df.apply(calculate_likelihood_absent, axis=1, args=(error_rate,))
    df['log_likelihood_present'] = df.apply(calculate_likelihood_present, axis=1, args=(error_rate,))

    return df


def calculate_likelihood_absent(row, e_s):
    return log_likelihood_absent(
        e_s,
        row['alt_counts'],
        row['alt_counts'] + row['ref_counts'],
    )


def calculate_likelihood_present(row, e_s):
    return log_likelihood_present(
        row['major_cn'],
        row['major_cn'] + row['minor_cn'],
        e_s,
        row['alt_counts'],
        row['ref_counts'] + row['alt_counts'],
    )

def log_binomial_pdf(x, n, p):
    return stats.binom.logpmf(x, n, p)


def log_likelihood_absent(e_s, n_v, n_t):
    return log_binomial_pdf(n_v, n_t, e_s)


def log_likelihood_present(c_m, c_t, e_s, n_v, n_t):
    if c_m == 0:
        return log_likelihood_absent(e_s, n_v, n_t)

    conditional_log_likelihoods = []

    for c_v in np.arange(1., c_m + 1., 1.):
        r = c_v / c_t
        conditional_log_likelihoods.append(log_binomial_pdf(n_v, n_t, r))

    return log_sum_exp(conditional_log_likelihoods)
