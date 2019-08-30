import numpy as np
from scgenome.constants import MAX_CN


def cn_mat_poisson(num_sample, num_bin, init_rng=np.random.poisson,
                   jump_rng=np.random.poisson, init_lambda=1., jump_lambda=1.,
                   seed=None, max_cn=MAX_CN):
    if seed is not None:
        np.random.seed(seed)
    first = init_rng(lam=init_lambda, size=num_sample)

    num_jump = (num_bin - 1) * num_sample
    if seed is not None:
        np.random.seed(seed)
    all_jumps = jump_rng(lam=jump_lambda, size=num_jump)

    if seed is not None:
        np.random.seed(seed)
    signs = np.random.binomial(n=1, p=0.5, size=num_jump)
    signs[np.where(signs == 0)] = -1

    cn_mat = np.zeros((num_sample, num_bin))
    cn_mat[:, 0] = first

    i = 0
    for r in range(0, num_sample):
        for c in range(1, num_bin):
            cn_mat[r, c] = max(cn_mat[r, c-1] + signs[i]*all_jumps[i], 0)
            cn_mat[r, c] = min(cn_mat[r, c], max_cn)
            i += 1

    return cn_mat.astype("int")


def get_poisson_bicluster(samples_per_cluster, num_bin, max_cn,
                          init_lambdas=(None, None),
                          jump_lambdas=(None, None), seeds=(None, None)):
    cluster1 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[0],
                              jump_lambda=jump_lambdas[0], seed=seeds[0],
                              max_cn=max_cn)
    cluster2 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[1],
                              jump_lambda=jump_lambdas[1], seed=seeds[1],
                              max_cn=max_cn)
    cn_mat = np.concatenate([cluster1, cluster2])
    return cn_mat