import numpy as np
import pandas as pd
from .constants import MAX_CN, INIT_CN, CHR_CLUSTERING, TRANS_NOT_SYMMETRIC, \
    TRANS_NOT_SUM, TRANS_NOT_SQUARE, NBIN_NCHR
from .utils import cn_mat_as_df

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

# TODO INIT_CN is 2
def clustered_cn_mat(num_sample, num_bin, chr_clustering, trans_mats,
                     max_cn=MAX_CN, init_cn=INIT_CN, seeds=None,
                     output_df=False, split=False):
    """

    :param num_sample: Number of samples/cells
    :param num_bin: How many bins the genome is separated into, must be
    divisible by number of chromosomes
    :param chr_clustering: array of ints describing clustering scheme for
    chromosomes. Eg. [0, 0, 1] means chromosomes 1 and 2 are in the same
    cluster, and chromosome 3 is in a different cluster. Elements of array
    index trans_mats.
    :param trans_mats: list of transition probability matrices to simulate
    copy-number profile with
    :return:
    """
    # TODO default trans mat
    if num_bin % len(chr_clustering) != 0:
        raise ValueError(NBIN_NCHR)
    if max(chr_clustering) > len(trans_mats) - 1:
        raise ValueError(CHR_CLUSTERING)

    clustered_trans_mats = [trans_mats[i] for i in chr_clustering]

    bins_per_chr = int(num_bin / len(chr_clustering))

    chr_mats = [cn_mat_from_trans(num_sample, bins_per_chr, trans, max_cn,
                                  init_cn, seeds)
                for trans in clustered_trans_mats]

    if split:
        return chr_mats

    cn_mat = np.concatenate(chr_mats, axis=1)
    if not output_df:
        return cn_mat
    else:
        chr_names = [str(c) for c in range(len(chr_clustering))]
        return cn_mat_as_df(cn_mat, chr_names)


def cn_mat_from_trans(num_sample, num_bin, trans, max_cn=MAX_CN,
                      init_cn=INIT_CN, seeds=None):
    cn_mat = np.zeros((num_sample, num_bin))
    for i in range(num_sample):
        cn_mat[i, :] = cell_from_trans(num_bin, trans, max_cn, init_cn, seeds)
    return cn_mat.astype("int")


def cell_from_trans(num_bin, trans, max_cn=MAX_CN, init_cn=INIT_CN,
                    seeds=None):
    trans_valid(trans)

    cn = np.zeros(num_bin)
    init_cn = min(init_cn, max_cn)
    cn[0] = init_cn

    for i in range(1, num_bin):
        choices = list(range(max_cn+1))
        probabilities = trans[cn[i-1].astype("int"), :]

        if seeds is not None:
            np.random.seed(seeds[i % len(seeds)])

        next_cn = np.random.choice(a=choices, size=1, p=probabilities)[0]
        cn[i] = next_cn

    return cn


def simulate_hmmcopy_data(num_sample, num_bin, num_chr, distribution="poisson",
                          **kwargs):
    if num_bin % num_chr != 0:
        raise ValueError("num_bin must be divisible by num_chr")

    if distribution == "poisson":
        cn_mat = cn_mat_poisson(num_sample, num_bin, **kwargs)
        # TODO create wrapper method so there is only one simulate_mat()

    num_bins_in_chr = int(num_bin / num_chr)

    m_cn_mat = pd.melt(pd.DataFrame(cn_mat), var_name="bin", value_name="copy")
    m_cn_mat.index.name = 'cell'
    # TODO double check that cell is cell and bin is bin
    m_cn_mat.reset_index(inplace=True)

    num_rep_chr = num_sample * num_bins_in_chr
    chr_vec = np.repeat(np.array(list(range(1, num_chr+1))), num_rep_chr)

    m_cn_mat["chr"] = chr_vec

    return m_cn_mat

def trans_valid(mat):
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(TRANS_NOT_SQUARE)
    if not check_symmetric(mat):
        raise ValueError(TRANS_NOT_SYMMETRIC)
    if not check_sum_valid(mat):
        raise ValueError(TRANS_NOT_SUM)

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)

def check_sum_valid(a):
    sum0 = np.sum(a, 0) == 1
    sum1 = np.sum(a, 1) == 1
    return np.all(sum0) and np.all(sum1)

def cn_mat_as_list(mat, chr_names):
    """
    Converts a sample x bin CN profile matrix to a pandas DataFrame to be
    plotted by scgenome.cnplot.plot_cell_cn_profile function in scgenome repo
    """
    column_names = ["chr", "start", "end", "copy1", "copy2"]
    df = cn_mat_as_df(mat, chr_names)
    melted = pd.melt(df, var_name="chr_bin", value_name="copy")
    chr_bin = melted["chr_bin"].str.split("_", expand=True)
    melted["chr"] = chr_bin[0]
    melted["bin"] = chr_bin[1]

    melted["cell_id"] = np.repeat(list(range(mat.shape[0])), mat.shape[1])
    melted = melted.drop("chr_bin", axis=1)
    melted = melted[["chr", "bin", "cell_id", "copy"]]
    return melted


