# Do BHC on a spike in experiment
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot
import scipy.cluster.hierarchy as sch
import os
import time
import matplotlib.pyplot as plt
import sklearn.metrics as skm
import sys

from scgenome.constants import LOG_P5

OUT_DIR = "/Users/massoudmaher/data/ntest_bhc/"
CN_DATA_FP = "/Users/massoudmaher/data/sc1935_clean_qc.csv"
#SAMPLE_IDS = ['SC-1935', 'SC-1936', 'SC-1937']
SAMPLE_IDS = ['SC-1935']
N_CELLS = 50  # Set to None if we want all cells
N_BIN = 100
spike_in = False
PROPORTIONS = None  # Set to None for equal proportion of each sample
N_STATES = 12
ALPHA = 10
PROB_CN_CHANGE = 0.8
params = pd.DataFrame({
    "CN_DATA_FP": CN_DATA_FP,
    "SAMPLE_IDS": SAMPLE_IDS,
    "N_CELLS": N_CELLS,
    "spike_in": spike_in,
    "PROPORTIONS": PROPORTIONS,
    "N_STATES": N_STATES,
    "ALPHA": ALPHA,
    "PROB_CN_CHANGE": PROB_CN_CHANGE
})

if len(sys.argv) >= 3:
    OUT_DIR = sys.argv[1]
    CN_DATA_FP = sys.argv[2]
    if len(sys.argv) >= 4:
        N_CELLS = int(sys.argv[3])
    if len(sys.argv) >= 5:
        N_BIN = int(sys.argv[4])

    print(f"args {OUT_DIR}, {CN_DATA_FP}, {N_CELLS} {N_BIN}")

if not os.path.exists(OUT_DIR):
    print(f"{OUT_DIR} does not exist, creating it")
    os.makedirs(OUT_DIR)
params.to_csv(os.path.join(OUT_DIR, "params.csv"))

print(f"Reading in {CN_DATA_FP}")
cn_data = pd.read_csv(CN_DATA_FP)

if N_CELLS is not None:
    if spike_in:
        cn_data = utils.get_cn_data_submixture(cn_data, N_CELLS, SAMPLE_IDS,
                                               proportions=PROPORTIONS)
    else:
        cn_data = utils.subsample_cn_data(cn_data, N_CELLS)

if N_BIN is not None:
    end_val = np.unique(cn_data["end"])[N_BIN]
    print(f"Reducing to {N_BIN} bins, end: {end_val}")
    cn_data = cn_data[cn_data["end"] <= end_val]

n_cell = np.unique(cn_data["cell_id"]).shape[0]
n_bin = np.unique(cn_data["end"]).shape[0]

print(f"Doing BHC on {n_cell} cells, {n_bin} bins")
start = time.time()
bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=N_STATES, alpha=ALPHA,
                               prob_cn_change=PROB_CN_CHANGE,
                               debug=True, clustering_id="copy")
)
print(f"{time.time()-start}s for BHC on {n_cell} cells, {n_bin} bins")

print(f"measurement {measurement.shape}, variances {variances.shape}")

bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)
lbhc_plot_data = bhc_plot_data.copy()
lbhc_plot_data[:, 2] = np.log(bhc_plot_data[:, 2])

bhc_linkage.to_csv(os.path.join(OUT_DIR, "linkage.csv"))
bhc_cell_ids.to_csv(os.path.join(OUT_DIR, "cell_ids.csv"), header=False)
matrix_data.to_csv(os.path.join(OUT_DIR, "matrix_data.csv"))
np.savetxt(os.path.join(OUT_DIR, "measurement.txt"), measurement)
np.savetxt(os.path.join(OUT_DIR, "variances.txt"), variances)
np.savetxt(os.path.join(OUT_DIR, "bhc_plot_data.txt"), bhc_plot_data)
np.savetxt(os.path.join(OUT_DIR, "lbhc_plot_data.txt"), lbhc_plot_data)

################ Plotting

#bhc_clusters = sch.fcluster(lbhc_plot_data, 12, criterion="distance")
bhc_clusters = sch.cut_tree(bhc_plot_data, height=-1*LOG_P5).flatten()
assert len(set(bhc_clusters)) > 1
cn_data = cncluster.prune_cluster(bhc_clusters, bhc_cell_ids, cn_data)

fig = plt.figure(figsize=(10, 8))
origin = "origin_id_int" if spike_in else None
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="bhc_cluster_id",
    linkage=lbhc_plot_data, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=bhc_cell_ids)
fig.savefig(os.path.join(OUT_DIR, "heatmap.pdf"), bbox_inches='tight')

################## Metrics
clabels = utils.get_mixture_labels(cn_data)
scores = skm.homogeneity_completeness_v_measure(clabels["origin_id_int"],
                                                clabels["bhc_cluster_id"])
print(f"homogeneity: {scores[0]}, completeness: {scores[1]}, "
      f"v-measure: {scores[2]}")

