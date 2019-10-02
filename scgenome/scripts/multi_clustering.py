# Do BHC, Naive and UMAP/HDBSCAN multiple times
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot, qc
import scipy.cluster.hierarchy as sch
from os.path import join
import os
import time
import sklearn.metrics as skm
import seaborn as sns

#OUT_DIR = "/Users/massoudmaher/data/test_mult_clst/"
#METRICS_FP = "/Users/massoudmaher/data/high_q_low_r/sc1661_1802_metrics.csv"
OUT_DIR = "/work/shah/maherm/low_reads100/"
METRICS_FP = "/work/shah/maherm/high_q_low_r/sc1661_1802_metrics.csv"
SEED = None

N_CELLS = 100  # Set to None if we want all cells
N_BIN = None
N_TRIAL = 30

# QC params
#N_READ = [50e3, 100e3, 250e3, 500e3]
N_READ = [300e3, 400e3, 500e3]
LIMIT_FPS = [
    #"sc1661_1802_qc_cnd_reads_lt_200k.csv",
    "sc1661_1802_qc_cnd_reads_lt_300k.csv",
    "sc1661_1802_qc_cnd_reads_lt_400k.csv",
    "sc1661_1802_qc_cnd_reads_lt_500k.csv"
]
#LIMIT_FPS = [join("/Users/massoudmaher/data/high_q_low_r", fp)
LIMIT_FPS = [join("/work/shah/maherm/high_q_low_r", fp)
             for fp in LIMIT_FPS]

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
#SAMPLE_IDS = ['SC-1935', 'SC-1936', 'SC-1937']
SAMPLE_IDS = ['SC-1661', 'SC-1802']
spike_in = True
PROPORTIONS = None  # Set to None for equal proportion of each sample

if not os.path.exists(OUT_DIR):
    print(f"{OUT_DIR} does not exist, creating it")
    os.makedirs(OUT_DIR)

if not os.path.exists(join(OUT_DIR, "plots")):
    os.makedirs(join(OUT_DIR, "plots"))

limit_cnds = []
for fp in LIMIT_FPS:
    cn_data = pd.read_csv(fp)
    cn_data["origin_id"] = cn_data["sample_id"]
    cn_data["origin_id_int"] = cn_data["origin_id"].factorize()[0]
    limit_cnds.append(cn_data)

#cnds = [qc.filt_reads(metrics, cn_data, n_reads=nr) for nr in N_READ]
cnds = []
N_READ = N_TRIAL*N_READ
limit_read_ind = list(range(len(limit_cnds))) * N_TRIAL
#for nr in N_READ:
for i in range(len(N_READ)):
    nr = N_READ[i]
    print(f"limit_read_int {limit_read_ind[i]}")
    print(f"fp[limit_read_int] {LIMIT_FPS[limit_read_ind[i]]}")
    print(f"nr {nr}")
    #n_cell = len(cnd["cell_id"].unique())
    cnd = utils.get_cn_data_submixture(limit_cnds[limit_read_ind[i]],
                                       N_CELLS, SAMPLE_IDS,
                                       proportions=PROPORTIONS)
    print(f"post cnd.shape {cnd.shape}")
    if N_BIN is not None:
        end_val = np.unique(cnd["end"])[N_BIN]
        print(f"Reducing to {N_BIN} bins, end: {end_val}")
        cnd = cnd[cnd["end"] <= end_val]

    if N_CELLS is not None or N_BIN is not None:
        cnd.to_csv(os.path.join(OUT_DIR, "plots", f"{nr}r_cn_data.csv"))
    cnds.append(cnd)


n_cell = N_CELLS
n_bin = -1

out_rows = []
for i in range(len(N_READ)):
    nr = N_READ[i]
    cnd = cnds[i]
    n_cell = len(set(cnd["cell_id"]))
    n_bin = len(set(cnd["end"]))
    print(f"n_read: {nr}")

    ######### BHC
    print(f"--Doing BHC on {n_cell} cells, {n_bin} bins")
    start = time.time()
    bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
        cncluster.bayesian_cluster(cnd, n_states=N_STATES, alpha=ALPHA,
                                   prob_cn_change=PROB_CN_SAME,
                                   debug=True, clustering_id="copy")
    )
    print(f"--{time.time()-start}s for BHC. {n_cell} cells, {n_bin} bins\n\n")
    bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)
    lbhc_plot_data = bhc_plot_data.copy()
    lbhc_plot_data[:, 2] = np.log(lbhc_plot_data[:, 2])

    #bhc_linkage.to_csv(join(OUT_DIR, f"{nr}reads_bhc_linkage.csv"))
    #bhc_cell_ids.to_csv(join(OUT_DIR, f"{nr}reads_bhc_cell_ids.csv"),
    #                    header=False)
    #np.savetxt(join(OUT_DIR, f"{nr}reads_bhc_plot_data.txt"), bhc_plot_data)

    # Try different thresholds
    grid = simulation.expand_grid({
        "transform": ["log"],
        "criterion": ["distance"],
        "threshold": np.arange(3, 15, step=1)
    })
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
    grid.to_csv(join(OUT_DIR, "plots", f"{nr}r_bhc_grid.csv"))

    # Plot grid for reference
    fig = plt.figure(figsize=(8, 8))
    g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                      sharey=True, sharex=False).add_legend()
    g = g.map(plt.scatter, "threshold", "num_clusters")
    plt.savefig(join(OUT_DIR, "plots", f"{nr}r_bhc_cluster_thresh.png"))
    plt.close()

    # Select threshold that maximizes metric
    i = grid["v_measure"].idxmax()
    bhc_max_row = grid.iloc[i, :]
    v_meas = grid["v_measure"][i]
    rand = grid["rand"][i]
    num_clusters = grid["num_clusters"][i]
    transform = grid["transform"][i]
    criterion = grid["criterion"][i]
    thresh = grid["threshold"][i]
    out_rows.append([nr, v_meas, rand, num_clusters,
                     "BHC", transform, criterion, thresh])

    ############### NHC
    print(f"--Doing NHC on {n_cell} cells, {n_bin} bins")
    start = time.time()
    naive_linkage = sch.linkage(np.nan_to_num(measurement),
                                method=NAIVE_METHOD, metric=NAIVE_METRIC)
    print(f"--{time.time()-start}s for NHC. {n_cell} cells, {n_bin} bins\n\n")
    lnaive_linkage = naive_linkage.copy()
    lnaive_linkage[:, 2] = np.log(lnaive_linkage[:, 2])
    #np.savetxt(join(OUT_DIR, "naive_plot_data.txt"), naive_linkage)

    # Try different thresholds
    grid = simulation.expand_grid({
        "transform": ["log"],
        "criterion": ["distance"],
        "threshold": np.arange(0, 5, step=0.25)
    })
    def apply_fn(row):
        if row["transform"] == "log":
            z = lnaive_linkage
        else:
            z = naive_linkage
        return sch.fcluster(z, row["threshold"], criterion=row["criterion"])
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
    grid.to_csv(join(OUT_DIR, "plots", f"{nr}r_naive_grid.csv"))

    # Plot grid for reference
    fig = plt.figure(figsize=(8, 8))
    g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                      sharey=True, sharex=False).add_legend()
    g = g.map(plt.scatter, "threshold", "num_clusters")
    plt.savefig(join(OUT_DIR, "plots", f"{nr}r_naive_cluster_thresh.png"))
    plt.close()

    # Select threshold that maximizes metric
    i = grid["v_measure"].idxmax()
    naive_max_row = grid.iloc[i, :]
    v_meas = grid["v_measure"][i]
    rand = grid["rand"][i]
    num_clusters = grid["num_clusters"][i]
    transform = grid["transform"][i]
    criterion = grid["criterion"][i]
    thresh = grid["threshold"][i]
    out_rows.append([nr, v_meas, rand, num_clusters,
                     "naive", transform, criterion, thresh])

    ############# UMAP/HDBSCAN
    print(f"--Doing UM+HDB on {n_cell} cells, {n_bin} bins")
    start = time.time()
    cn = (cnd.set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0))
    uh_cluster = cncluster.umap_hdbscan_cluster(cn, n_components=2,
                                                n_neighbors=UMAP_NN,
                                                min_dist=UMAP_MIN_DIST)
    print(f"--{time.time()-start}s for UMHD. {n_cell} cells, {n_bin} bins\n\n")
    uh_cluster.columns = ['cell_id', 'umap_cluster_id', 'umap1', 'umap2']
    uh_cluster.to_csv(join(OUT_DIR, "plots", f"{nr}r_umap_hdb_clust.csv"))
    cnd = cnd.merge(uh_cluster)

    v_meas = skm.homogeneity_completeness_v_measure(cnd["origin_id_int"],
                                                    cnd["umap_cluster_id"])[2]
    rand = skm.adjusted_rand_score(cnd["origin_id_int"],
                                   cnd["umap_cluster_id"])
    out_rows.append([nr, v_meas, rand, None, "umap", None, None, None])

    out_mat = pd.DataFrame(out_rows,
        columns=["max_n_reads", "v_measure", "rand_index", "num_clusters",
                 "method", "transform", "criterion", "threshold"])
    out_mat.to_csv(join(OUT_DIR, "out_mat.csv"))

fig = plt.figure(figsize=(8, 8))
sns.boxplot(data=out_mat, x="max_n_reads", y="v_measure", hue="method")
plt.savefig(join(OUT_DIR, "v_measure.png"))

fig = plt.figure(figsize=(8, 8))
sns.boxplot(data=out_mat, x="max_n_reads", y="rand_index", hue="method")
plt.savefig(join(OUT_DIR, "rand.png"))
