import numpy as np
import pandas as pd
from tqdm import tqdm
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from scipy.cluster.hierarchy import linkage

from scgenome import cncluster
from scgenome.constants import SIM_META, REQUIRED_PARAMS, NAIVE_METRIC
from scgenome.simulation.gaussian import get_gaussian_bicluster
from scgenome.simulation.poisson import get_poisson_bicluster
from scgenome.utils import cn_mat_as_df, cn_mat_to_cn_data, expand_grid


def do_simulations(params, naive_metric=NAIVE_METRIC, num_cores=None):
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
        simulate_set(samples_per_cluster, num_bin, max_cn, alpha,
                     distribution, dist_params, df)
        cn_mat = df["cn_mat"]
        cn_data = df["cn_data"]
        cell_ids = df["cell_ids"]
        do_bhc(cn_mat, cn_data, max_cn, alpha, cell_ids, df=df)
        do_naive_hc(cn_mat, metric=naive_metric, df=df)

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
        cn_mat = get_poisson_bicluster(samples_per_cluster, num_bin, max_cn,
                                       df=None, init_lambdas=init_lambdas,
                                       jump_lambdas=jump_lambdas, seeds=seeds,
                                       noise_seed=noise_seed)
    elif distribution == "gaussian":
        cn_mat = get_gaussian_bicluster(samples_per_cluster, num_bin, max_cn)
    else:
        raise ValueError(f"distribution {distribution} invalid")
    _bi_cn_mat_assignment(cn_mat, samples_per_cluster, df)


def _bi_cn_mat_assignment(cn_mat, samples_per_cluster, df=None):
    clst1_cell_ids = [f"cl1_cell{i}" for i in range(samples_per_cluster)]
    clst2_cell_ids = [f"cl2_cell{i}" for i in range(samples_per_cluster)]
    cell_ids = clst1_cell_ids + clst2_cell_ids

    chr_names = ["1", "2"]
    df_cn_mat = cn_mat_as_df(cn_mat, chr_names)
    cn_data = cn_mat_to_cn_data(df_cn_mat, cell_id_vals=cell_ids)
    cn_data["cluster_id"] = (
        cn_data["cell_id"].str.split("_", expand=True).iloc[:, 0])
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


def do_naive_hc(cn_mat, metric=NAIVE_METRIC, df=None):
    link_mat = linkage(cn_mat, metric=metric)
    if df is not None:
        df["naive_hc"] = link_mat
        return df
    else:
        return link_mat


def get_prop_correct(clustering):
    return max((clustering["exp_cl"] == clustering["obs_cl"]).value_counts() /
               clustering.shape[0])


def get_plot_data(linkage):
    plinkage = linkage.loc[:, ["i", "j", "r_merge", "merge_count"]]
    plinkage["r_merge"] = plinkage["r_merge"].astype("float")
    plinkage["dist"] = -1 * plinkage["r_merge"]
    plot_data = (
        plinkage[["i", "j", "dist", "merge_count"]].to_numpy().astype("float"))
    return plinkage, plot_data


