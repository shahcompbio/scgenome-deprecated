# Do BHC, Naive and UMAP/HDBSCAN multiple times
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot, qc
import scipy.cluster.hierarchy as sch
from os.path import join
import os
import time
import sklearn.metrics as skm
import matplotlib.pyplot as plt
import seaborn as sns
import random

OUT_DIR = "/Users/massoudmaher/data/test_vary_prop2/"
CN_DATA_FP = "/Users/massoudmaher/data/clean_sc_1935_1936_1937_cn_data_qc.csv"
SEED = None

N_CELLS = 100  # Set to None if we want all cells
N_BIN = None

# BHC Params
N_STATES = 12
ALPHA = 10
PROB_CN_SAME = 0.8

# Naive params
NAIVE_METHOD = 'single'
NAIVE_METRIC = 'euclidean'

# UMAP / HDBSCAN params
UMAP_NN = 5
UMAP_MIN_DIST = 0.1

# Spike in params (not always used)
SAMPLE_IDS = ['SC-1935', 'SC-1936', 'SC-1937']
spike_in = True
PROPORTIONS = [
    [0.05, 0.35, 0.6],
    [0.1, 0.3, 0.6],
    [0.2, 0.4, 0.4]
    #[0.3, 0.4, 0.4],
]  # Set to None for equal proportion of each sample
TRIAL_PER_PROP = 20

if not os.path.exists(OUT_DIR):
    print(f"{OUT_DIR} does not exist, creating it")
    os.makedirs(OUT_DIR)

if not os.path.exists(join(OUT_DIR, "plots")):
    os.makedirs(join(OUT_DIR, "plots"))

print(f"Reading in {CN_DATA_FP}")
cn_data = pd.read_csv(CN_DATA_FP)

if N_BIN is not None:
    end_val = np.unique(cn_data["end"])[N_BIN]
    print(f"Reducing to {N_BIN} bins, end: {end_val}")
    cn_data = cn_data[cn_data["end"] <= end_val]

prop_prod = PROPORTIONS*TRIAL_PER_PROP
print(f"prop_prod {prop_prod}")
cnds = []
for prop in prop_prod:
    cnd = utils.get_cn_data_submixture(cn_data, N_CELLS, SAMPLE_IDS,
                                       proportions=prop)
    random.shuffle(SAMPLE_IDS)
    cnds.append(cnd)


out_rows = []
for i in range(len(cnds)):
    cnd = cnds[i]
    n_cell = len(set(cnd["cell_id"]))
    n_bin = len(set(cnd["end"]))
    prop = prop_prod[i]
    print(f"prop: {prop}")

    ############ BHC
    print(f"--Doing BHC on {n_cell} cells, {n_bin} bins")
    start = time.time()
    bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
        cncluster.bayesian_cluster(cnd, n_states=N_STATES, alpha=ALPHA,
                                   prob_cn_change=PROB_CN_SAME,
                                   debug=True, clustering_id="copy")
    )
    print(f"--{time.time()-start}s for BHC. {n_cell} cells, {n_bin} bins\n\n")
    bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)

    grid = simulation.expand_grid({
        "transform": ["log"],
        "criterion": ["distance"],
        "threshold": np.arange(0, 20, step=1)
    })

    lbhc_plot_data = bhc_plot_data.copy()
    lbhc_plot_data[:, 2] = np.log(lbhc_plot_data[:, 2])

    # Try different thresholds
    def apply_fn(row):
        if row["transform"] == "log":
            z = lbhc_plot_data
        else:
            z = bhc_plot_data
        return sch.fcluster(z, row["threshold"], criterion=row["criterion"])
    grid["fcluster"] = grid.apply(apply_fn, axis=1)
    grid["num_clusters"] = grid["fcluster"].apply(lambda x: len(set(x)))

    # Get metrics for each threshold
    cluster_id = "bhc_cluster_id"
    def apply_fn(row):
        clst_cnd = cncluster.prune_cluster(row["fcluster"], bhc_cell_ids, cnd,
                                           cluster_id)
        clabels = utils.get_mixture_labels(clst_cnd, obs_name=cluster_id)
        hcv = skm.homogeneity_completeness_v_measure(clabels["origin_id_int"],
                                                     clabels[cluster_id])
        rand = skm.adjusted_rand_score(clabels["origin_id_int"],
                                       clabels[cluster_id])
        return list(hcv) + [rand]
    scores = grid.apply(apply_fn, axis=1).tolist()
    scores = pd.DataFrame(scores, columns=["homogeneity", "completeness", "v_measure", "rand"])
    grid = pd.concat([grid, scores], axis=1)

    grid.to_csv(join(OUT_DIR, "plots", f"{i}i_bhc_grid.csv"))

    # Plot grid for reference
    fig = plt.figure(figsize=(8, 8))
    g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                      sharey=True, sharex=False).add_legend()
    g = g.map(plt.scatter, "threshold", "num_clusters")
    plt.savefig(join(OUT_DIR, "plots", f"{i}i_bhc_cluster_thresh.png"))

    # Select threshold that maximizes metric
    ind = grid["v_measure"].idxmax()
    bhc_max_row = grid.iloc[ind, :]
    v_meas = grid["v_measure"][ind]
    rand = grid["rand"][ind]
    min_prop = min(prop)
    out_rows.append([min_prop, v_meas, rand, "BHC"])

    ############# UMAP
    print(f"--Doing UM+HDB on {n_cell} cells, {n_bin} bins")
    start = time.time()
    cn = (cnd.set_index(['chr', 'start', 'end', 'cell_id'])['copy']
          .unstack(level='cell_id').fillna(0))
    uh_cluster = cncluster.umap_hdbscan_cluster(cn, n_components=2,
                                                n_neighbors=UMAP_NN,
                                                min_dist=UMAP_MIN_DIST)
    print(f"--{time.time()-start}s for UMHD. {n_cell} cells, {n_bin} bins\n\n")
    uh_cluster.columns = ['cell_id', 'umap_cluster_id', 'umap1', 'umap2']
    uh_cluster.to_csv(join(OUT_DIR, "plots", f"{i}i_umap_hdb_clust.csv"))
    cnd = cnd.merge(uh_cluster)
    v_meas = skm.homogeneity_completeness_v_measure(cnd["origin_id_int"],
                                                    cnd["umap_cluster_id"])[2]
    rand = skm.adjusted_rand_score(cnd["origin_id_int"],
                                   cnd["umap_cluster_id"])
    out_rows.append([min_prop, v_meas, rand, "umap"])

    ########### NHC
    print(f"--Doing NHC on {n_cell} cells, {n_bin} bins")
    start = time.time()
    naive_linkage = sch.linkage(np.nan_to_num(measurement),
                                method=NAIVE_METHOD, metric=NAIVE_METRIC)
    print(f"--{time.time()-start}s for NHC. {n_cell} cells, {n_bin} bins\n\n")
    #np.savetxt(join(OUT_DIR, "naive_plot_data.txt"), naive_linkage)

    # Try different thresholds
    grid = simulation.expand_grid({
        "transform": ["none"],
        "criterion": ["distance"],
        "threshold": np.arange(0, 20, step=1)
    })
    def apply_fn(row):
        return sch.fcluster(naive_linkage, row["threshold"],
                            criterion=row["criterion"])
    grid["fcluster"] = grid.apply(apply_fn, axis=1)
    grid["num_clusters"] = grid["fcluster"].apply(lambda x: len(set(x)))

    # Get metrics for each threshold
    cluster_id = "naive_cluster_id"
    def apply_fn(row):
        clst_cnd = cncluster.prune_cluster(row["fcluster"], bhc_cell_ids, cnd,
                                           cluster_id)
        clabels = utils.get_mixture_labels(clst_cnd, obs_name=cluster_id)
        hcv = skm.homogeneity_completeness_v_measure(clabels["origin_id_int"],
                                                     clabels[cluster_id])
        rand = skm.adjusted_rand_score(clabels["origin_id_int"],
                                       clabels[cluster_id])
        return list(hcv) + [rand]
    scores = grid.apply(apply_fn, axis=1).tolist()
    scores = pd.DataFrame(scores, columns=["homogeneity", "completeness", "v_measure", "rand"])
    grid = pd.concat([grid, scores], axis=1)
    grid.to_csv(join(OUT_DIR, "plots", f"{i}i_naive_grid.csv"))
    # Plot grid for reference
    fig = plt.figure(figsize=(8, 8))
    g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                      sharey=True, sharex=False).add_legend()
    g = g.map(plt.scatter, "threshold", "num_clusters")
    plt.savefig(join(OUT_DIR, "plots", f"{i}i_naive_cluster_thresh.png"))

    # Select threshold that maximizes metric
    ind = grid["v_measure"].idxmax()
    naive_max_row = grid.iloc[ind, :]
    v_meas = grid["v_measure"][ind]
    rand = grid["rand"][ind]
    min_prop = min(prop)
    out_rows.append([min_prop, v_meas, rand, "naive"])



    out_mat = pd.DataFrame(out_rows,
                           columns=["min_prop", "v_measure", "rand_index", "method"])
    out_mat.to_csv(join(OUT_DIR, "out_mat.csv"))





