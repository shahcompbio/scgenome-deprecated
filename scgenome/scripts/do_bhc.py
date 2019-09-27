# Do BHC on a spike in experiment
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot
import scipy.cluster.hierarchy as sch
import os
import time
import matplotlib.pyplot as plt
import sklearn.metrics as skm


OUT_DIR = "/Users/massoudmaher/data/test_do_bhc/"
CN_DATA_FP = "/Users/massoudmaher/data/clean_sc_1935_1936_1937_cn_data_qc.csv"
SAMPLE_IDS = ['SC-1935', 'SC-1936', 'SC-1937']
N_CELLS = 20  # Set to None if we want all cells
PROPORTIONS = None  # Set to None for equal proportion of each sample
N_STATES = 12
ALPHA = 10
PROB_CN_CHANGE = 0.8

if not os.path.exists(CN_DATA_FP):
    print(f"{CN_DATA_FP} does not exist, creating it")
    os.makedirs(CN_DATA_FP)

print(f"Reading in {CN_DATA_FP}")
cn_data = pd.read_csv(CN_DATA_FP)

if N_CELLS is not None:
    cn_data = utils.get_cn_data_submixture(cn_data, N_CELLS, SAMPLE_IDS,
                                           proportions=PROPORTIONS)

print(f"Doing BHC on {N_CELLS} cells")

start = time.time()
bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=N_STATES, alpha=ALPHA,
                               prob_cn_change=PROB_CN_CHANGE)
)
print(f"{time.time()-start}s for BHC on {N_CELLS} cells")

bhc_linkage.to_csv(os.path.join(OUT_DIR, "linkage.csv"))
bhc_cell_ids.to_csv(os.path.join(OUT_DIR, "cell_ids.csv"), header=False)
matrix_data.to_csv(os.path.join(OUT_DIR, "matrix_data.csv"))
np.savetxt(os.path.join(OUT_DIR, "measurement.txt"), measurement)
np.savetxt(os.path.join(OUT_DIR, "variances.txt"), variances)

################ Plotting
bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)
lbhc_plot_data = bhc_plot_data.copy()
lbhc_plot_data[:, 2] = np.log(bhc_plot_data[:, 2])

bhc_clusters = sch.fcluster(lbhc_plot_data, 12, criterion="distance")
assert len(set(bhc_clusters)) > 1
cn_data = cncluster.prune_cluster(bhc_clusters, bhc_cell_ids, cn_data)

fig = plt.figure(figsize=(10, 8))
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="bhc_cluster_id",
    linkage=lbhc_plot_data, origin_field_name="origin_id_int", raw=True,
    flip=True, cell_id_order=bhc_cell_ids)
fig.savefig(os.path.join(OUT_DIR, "heatmap.pdf"), bbox_inches='tight')

################## Metrics
clabels = utils.get_mixture_labels(cn_data)
scores = skm.homogeneity_completeness_v_measure(clabels["origin_id_int"],
                                                clabels["bhc_cluster_id"])
print(f"homogeneity: {scores[0]}, completeness: {scores[1]}, "
      f"v-measure: {scores[2]}")

