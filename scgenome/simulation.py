import numpy as np
import pandas as pd
from tqdm import tqdm
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

from scgenome import cncluster
from .constants import MAX_CN, INIT_CN, CHR_CLUSTERING, \
    TRANS_NOT_SYMMETRIC, TRANS_NOT_SUM, TRANS_NOT_SQUARE, NBIN_NCHR, SIM_META, \
    REQUIRED_PARAMS
from .utils import cn_mat_as_df, cn_mat_to_cn_data, expand_grid


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


def get_prop_correct(clustering):
    return max((clustering["exp_cl"] == clustering["obs_cl"]).value_counts() /
               clustering.shape[0])


def poisson_bicluster(samples_per_cluster, num_bin, max_cn, alpha, df=None,
                      init_lambdas=(None, None),
                      jump_lambdas=(None, None), seeds=(None, None),
                      noise_seed=None):
    cluster1 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[0],
                              jump_lambda=jump_lambdas[0], seed=seeds[0],
                              max_cn=max_cn)
    cluster2 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[1],
                              jump_lambda=jump_lambdas[1], seed=seeds[1],
                              max_cn=max_cn)

    clst1_cell_ids = [f"cl1_cell{i}" for i in range(samples_per_cluster)]
    clst2_cell_ids = [f"cl2_cell{i}" for i in range(samples_per_cluster)]

    cn_mat = np.concatenate([cluster1, cluster2])
    cell_ids = clst1_cell_ids + clst2_cell_ids

    chr_names = ["1", "2"]
    df_cn_mat = cn_mat_as_df(cn_mat, chr_names)
    cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=cell_ids)
    cn_data["cluster_id"] = (
        cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
    if noise_seed is not None:
        np.random.seed(noise_seed)
    cn_data["copy2"] = cn_data["copy"] + np.absolute(
        np.random.normal(size=cn_data.shape[0], scale=0.3))
    cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                       "cluster_id", "copy"]

    tlinkage, root, cl_cell_ids = (
        cncluster.bayesian_cluster(cn_data, n_states=max_cn,
                                   value_ids=["copy"], alpha=alpha))

    plinkage, plot_data = get_plot_data(tlinkage)

    clustering = pd.DataFrame()
    clustering["sample_inds"] = list(range(cn_mat.shape[0]))
    clustering["cell_id"] = cell_ids
    clustering["exp_cl"] = clustering["cell_id"].str[2]

    left_samples = [x.sample_inds[0] for x in root.left_child.get_leaves()]
    right_samples = [x.sample_inds[0] for x in root.right_child.get_leaves()]

    def fn(ind):
        if ind in left_samples:
            return "1"
        elif ind in right_samples:
            return "2"

    clustering["obs_cl"] = clustering["sample_inds"].apply(fn)

    prop_correct = get_prop_correct(clustering)

    if df is not None:
        df["cn_data"] = cn_data
        df["plinkage"] = plinkage
        df["plot_data"] = plot_data
        df["clustering"] = clustering
        df["prop_correct"] = prop_correct
        return df
    else:
        return cn_data, plinkage, plot_data, clustering, prop_correct


def get_plot_data(linkage):
    plinkage = linkage.loc[:, ["i", "j", "r_merge", "merge_count"]]
    plinkage["r_merge"] = plinkage["r_merge"].astype("float")
    plinkage["dist"] = -1 * plinkage["r_merge"]
    plot_data = (
        plinkage[["i", "j", "dist", "merge_count"]].to_numpy().astype("float"))
    return plinkage, plot_data

def many_poisson_bicluster(trials_per_set, samples_per_cluster, num_bin,
                           max_cn, alpha, init_lambdas, jump_lambdas,
                           num_cores=None):
    params = {"samples_per_cluster": samples_per_cluster,
              "num_bin": num_bin, "max_cn": max_cn, "alpha": alpha,
              "init_lambdas": init_lambdas, "jump_lambdas": jump_lambdas}
    sim_df = expand_grid(params)
    sim_df = pd.concat([sim_df] * trials_per_set)

    def apply_fn(df):
        samples_per_cluster = df["samples_per_cluster"]
        num_bin = df["num_bin"]
        max_cn = df["max_cn"]
        alpha = df["alpha"]
        init_lambdas = df["init_lambdas"]
        jump_lambdas = df["jump_lambdas"]
        return poisson_bicluster(samples_per_cluster, num_bin, max_cn, alpha,
                                 df=df, init_lambdas=init_lambdas,
                                 jump_lambdas=jump_lambdas)

    if num_cores is None:
        tqdm.pandas()
        sim_df = sim_df.progress_apply(apply_fn, axis=1).reset_index(drop=True)
        return sim_df
    else:
        ProgressBar().register()
        sim_df = dd.from_pandas(sim_df, npartitions=num_cores)
        result = sim_df.map_partitions(lambda df: df.apply(apply_fn, axis=1),
                                       meta=SIM_META)
        return result.compute(scheduler="processes").reset_index(drop=True)


def do_simulations(params):
    """"
    :param param_grid: dictionary of {<argument_name>: <arg_values> which is
    used to make parameter grid for simulation
    """
    if not REQUIRED_PARAMS.issubset(params.keys()):
        raise ValueError(f"param_grid requires these keys: {REQUIRED_PARAMS}")

    sim_df = expand_grid(params)
    sim_df = pd.concat([sim_df] * params["trials_per_set"])

    def apply_fn(df):
        samples_per_cluster = df["samples_per_cluster"]
        num_bin = df["num_bin"]
        max_cn = df["max_cn"]
        alpha = df["alpha"]
        distribution = df["distribution"]
        dist_params = df["distribution_params"]
        return poisson_bicluster(samples_per_cluster, num_bin, max_cn, alpha,
                                 df=df, init_lambdas=init_lambdas,
                                 jump_lambdas=jump_lambdas)

    if num_cores is None:
        tqdm.pandas()
        sim_df = sim_df.progress_apply(apply_fn, axis=1).reset_index(drop=True)
        return sim_df
    else:
        ProgressBar().register()
        sim_df = dd.from_pandas(sim_df, npartitions=num_cores)
        result = sim_df.map_partitions(lambda df: df.apply(apply_fn, axis=1),
                                       meta=SIM_META)
        return result.compute(scheduler="processes").reset_index(drop=True)


def simulate_set(samples_per_cluster, num_bin, max_cn, alpha, distribution,
                 dist_params, df=None):
    if distribution == "poisson_rw":
        init_lambdas = dist_params["init_lambdas"]
        jump_lambdas = dist_params["init_lambdas"]
        try:
            seeds = dist_params["seeds"]
        except KeyError:
            seeds = None
        try:
            noise_seed = dist_params["noise_seed"]
        except KeyError:
            noise_seed = None
        get_poisson_bicluster(samples_per_cluster, num_bin, max_cn, df=df,
                              init_lambdas=init_lambdas,
                              jump_lambdas=jump_lambdas, seeds=seeds,
                              noise_seed=noise_seed)
    elif distribution == "gaussian":
        pass
    else:
        raise ValueError(f"distribution {distribution} invalid")


def get_poisson_bicluster(samples_per_cluster, num_bin, max_cn, df=None,
                          init_lambdas=(None, None),
                          jump_lambdas=(None, None), seeds=(None, None),
                          noise_seed=None):
    cluster1 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[0],
                              jump_lambda=jump_lambdas[0], seed=seeds[0],
                              max_cn=max_cn)
    cluster2 = cn_mat_poisson(samples_per_cluster, num_bin,
                              init_lambda=init_lambdas[1],
                              jump_lambda=jump_lambdas[1], seed=seeds[1],
                              max_cn=max_cn)

    clst1_cell_ids = [f"cl1_cell{i}" for i in range(samples_per_cluster)]
    clst2_cell_ids = [f"cl2_cell{i}" for i in range(samples_per_cluster)]

    cn_mat = np.concatenate([cluster1, cluster2])
    cell_ids = clst1_cell_ids + clst2_cell_ids

    chr_names = ["1", "2"]
    df_cn_mat = cn_mat_as_df(cn_mat, chr_names)
    cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=cell_ids)
    cn_data["cluster_id"] = (
        cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
    if noise_seed is not None:
        np.random.seed(noise_seed)
    cn_data["copy2"] = cn_data["copy"] + np.absolute(
        np.random.normal(size=cn_data.shape[0], scale=0.3))
    cn_data.columns = ["chr", "bin", "cell_id", "state", "start", "end",
                       "cluster_id", "copy"]
    if df is not None:
        df["cn_mat"] = cn_mat
        df["cn_data"] = cn_data
        df["cell_ids"] = cell_ids
        return df
    else:
        return cn_mat, cn_data, cell_ids


def do_bhc(cn_mat, cn_data, max_cn, alpha, cell_ids, df=None):
    tlinkage, root, cl_cell_ids = (
        cncluster.bayesian_cluster(cn_data, n_states=max_cn,
                                   value_ids=["copy"], alpha=alpha))

    plinkage, plot_data = get_plot_data(tlinkage)

    clustering = pd.DataFrame()
    clustering["sample_inds"] = list(range(cn_mat.shape[0]))
    clustering["cell_id"] = cell_ids
    clustering["exp_cl"] = clustering["cell_id"].str[2]

    left_samples = [x.sample_inds[0] for x in root.left_child.get_leaves()]
    right_samples = [x.sample_inds[0] for x in root.right_child.get_leaves()]

    def fn(ind):
        if ind in left_samples:
            return "1"
        elif ind in right_samples:
            return "2"

    clustering["obs_cl"] = clustering["sample_inds"].apply(fn)

    prop_correct = get_prop_correct(clustering)

    if df is not None:
        df["cn_data"] = cn_data
        df["plinkage"] = plinkage
        df["plot_data"] = plot_data
        df["clustering"] = clustering
        df["prop_correct"] = prop_correct
        return df
    else:
        return cn_data, plinkage, plot_data, clustering, prop_correct


def poisson_bicluster(samples_per_cluster, num_bin, max_cn, alpha, df=None,
                      init_lambdas=(None, None),
                      jump_lambdas=(None, None), seeds=(None, None),
                      noise_seed=None):
    get_poisson_bicluster(samples_per_cluster, num_bin, max_cn, df=df,
                          init_lambdas=init_lambdas,
                          jump_lambdas=jump_lambdas, seeds=seeds,
                          noise_seed=noise_seed)
    cn_mat = df["cn_mat"]
    cn_data = df["cn_data"]
    cell_ids = df["cell_ids"]

    do_bhc(cn_mat, cn_data, max_cn, alpha, cell_ids, df=df)



