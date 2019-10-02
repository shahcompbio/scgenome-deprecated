# Do BHC, Naive and UMAP/HDBSCAN clustering on a dataset
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot
import scipy.cluster.hierarchy as sch
import os
import time
import matplotlib.pyplot as plt
import sklearn.metrics as skm

from scgenome.constants import LOG_P5

OUT_DIR = "/Users/massoudmaher/data/test/"
CN_DATA_FP = "/Users/massoudmaher/data/clean_sc_1935_1936_1937_cn_data_qc.csv"
SEED = None

N_CELLS = 150  # Set to None if we want all cells
N_BIN = None

# BHC Params
N_STATES = 12
ALPHA = 10
PROB_CN_SAME = 0.8

# Naive params
NAIVE_METHOD = 'single'
NAIVE_METRIC = 'euclidean'

# UMAP / HDBSCAN params
UMAP_NN = 15
UMAP_MIN_DIST = 0.1

# Spike in params (not always used)
SAMPLE_IDS = ['SC-1935', 'SC-1936', 'SC-1937']
spike_in = True
PROPORTIONS = None  # Set to None for equal proportion of each sample

if not os.path.exists(OUT_DIR):
    print(f"{OUT_DIR} does not exist, creating it")
    os.makedirs(OUT_DIR)

print(f"Reading in {CN_DATA_FP}")
cn_data = pd.read_csv(CN_DATA_FP)

if N_CELLS is not None:
    if spike_in:
        cn_data = utils.get_cn_data_submixture(cn_data, N_CELLS, SAMPLE_IDS,
                                               proportions=PROPORTIONS)
    else:
        cells = pd.Series(np.unique(cn_data["cell_id"]))
        keep_cells = cells.sample(N_CELLS, random_state=SEED)
        cn_data = cn_data[cn_data["cell_id"].isin(keep_cells)]

if N_BIN is not None:
    end_val = np.unique(cn_data["end"])[N_BIN]
    print(f"Reducing to {N_BIN} bins, end: {end_val}")
    cn_data = cn_data[cn_data["end"] <= end_val]

if N_CELLS is not None or N_BIN is not None:
    cn_data.to_csv(os.path.join(OUT_DIR, "cn_data.csv"))

n_cell = np.unique(cn_data["cell_id"]).shape[0]
n_bin = np.unique(cn_data["end"]).shape[0]

print(f"cn_data.shape {cn_data.shape}")

print(f"Doing BHC on {n_cell} cells, {n_bin} bins")
start = time.time()
bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=N_STATES, alpha=ALPHA,
                               prob_cn_change=PROB_CN_SAME,
                               debug=True, clustering_id="copy")
)
print(f"{time.time()-start}s for BHC on {n_cell} cells, {n_bin} bins\n\n")
bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)
lbhc_plot_data = bhc_plot_data.copy()
lbhc_plot_data[:, 2] = np.log(bhc_plot_data[:, 2])
bhc_linkage.to_csv(os.path.join(OUT_DIR, "bhc_linkage.csv"))
bhc_cell_ids.to_csv(os.path.join(OUT_DIR, "bhc_cell_ids.csv"), header=False)
np.savetxt(os.path.join(OUT_DIR, "measurement.txt"), measurement)
np.savetxt(os.path.join(OUT_DIR, "variances.txt"), variances)
np.savetxt(os.path.join(OUT_DIR, "bhc_plot_data.txt"), bhc_plot_data)
np.savetxt(os.path.join(OUT_DIR, "lbhc_plot_data.txt"), lbhc_plot_data)

print(f"Doing NHC on {n_cell} cells, {n_bin} bins")
start = time.time()
naive_linkage = sch.linkage(np.nan_to_num(measurement), method=NAIVE_METHOD,
                            metric=NAIVE_METRIC)
print(f"{time.time()-start}s for NHC on {n_cell} cells, {n_bin} bins\n\n")
np.savetxt(os.path.join(OUT_DIR, "naive_plot_data.txt"), naive_linkage)

print(f"Doing UM+HDB on {n_cell} cells, {n_bin} bins")
start = time.time()
cn = (cn_data.set_index(['chr', 'start', 'end', 'cell_id'])['copy']
        .unstack(level='cell_id').fillna(0))
uh_cluster = cncluster.umap_hdbscan_cluster(cn, n_components=2,
                                            n_neighbors=UMAP_NN,
                                            min_dist=UMAP_MIN_DIST)
print(f"{time.time()-start}s for UM+HDB on {n_cell} cells, {n_bin} bins\n\n")
uh_cluster.columns = ['cell_id', 'umap_cluster_id', 'umap1', 'umap2']
uh_cluster.to_csv(os.path.join(OUT_DIR, "umap_hdb_clust.csv"))

print("Plotting")
# BHC
bhc_clusters = sch.fcluster(lbhc_plot_data, 12, criterion="distance")
assert len(set(bhc_clusters)) > 1
cn_data = cncluster.prune_cluster(bhc_clusters, bhc_cell_ids, cn_data,
                                  cluster_field_name="bhc_cluster_id")

origin = 'origin_id_int' if spike_in else None

fig = plt.figure(figsize=(10, 8))
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="bhc_cluster_id",
    linkage=lbhc_plot_data, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=bhc_cell_ids)
fig.savefig(os.path.join(OUT_DIR, "bhc_heatmap.png"), bbox_inches='tight')

# NHC
naive_clusters = sch.fcluster(naive_linkage, 12, criterion="distance")
assert len(set(naive_clusters)) > 1
cn_data = cncluster.prune_cluster(naive_clusters, bhc_cell_ids, cn_data,
                                  cluster_field_name="naive_cluster_id")

fig = plt.figure(figsize=(10, 8))
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="naive_cluster_id",
    linkage=naive_linkage, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=bhc_cell_ids)
fig.savefig(os.path.join(OUT_DIR, "naive_heatmap.png"), bbox_inches='tight')

# UMAP+HDBSCAN
cn_data = cn_data.merge(uh_cluster)
# Scatterplot
fig = plt.figure(figsize=(8, 8))
nuh_cluster = uh_cluster.copy()
nuh_cluster.columns = ['cell_id', 'cluster_id', 'umap1', 'umap2']
cncluster.plot_umap_clusters(plt.gca(), nuh_cluster)
fig.savefig(os.path.join(OUT_DIR, "uh_scatter.png"), bbox_inches='tight')
# Heatmap
fig = plt.figure(figsize=(10, 8))
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="umap_cluster_id",
    linkage=None, origin_field_name=origin, raw=True,
    flip=False)
fig.savefig(os.path.join(OUT_DIR, "umap_heatmap.png"), bbox_inches='tight')

# TODO save cn_data with all the clustering columns

################## Metrics
clabels = utils.get_mixture_labels(cn_data)
scores = skm.homogeneity_completeness_v_measure(clabels["origin_id_int"],
                                                clabels["bhc_cluster_id"])
print(f"homogeneity: {scores[0]}, completeness: {scores[1]}, "
      f"v-measure: {scores[2]}")

