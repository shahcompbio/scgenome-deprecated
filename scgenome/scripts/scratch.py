from scgenome import tantalus
import pandas as pd
from IPython.display import display
from scgenome import utils, cncluster, simulation, cnplot, jointcnmodels
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import sklearn.metrics as skm

def get_mixture_labels(cn_data, obs_name="bhc_cluster_id", 
                       exp_name="origin_id_int"):
    sub_cn_data = cn_data[["cell_id", obs_name, exp_name]].drop_duplicates()
    return sub_cn_data

cn_data_fp = "/Users/massoudmaher/data/sc1935to37_spike_in_chrX.csv"
cn_data = pd.read_csv(cn_data_fp).iloc[:,1:]

cn_data.head(10)
cn_data = cn_data[cn_data["end"] <= 5000000]

hmmcopy_tickets = ['SC-1935', 'SC-1936', 'SC-1937']
sample_ids = [["SA922"], ['SA921'], ['SA1090']]

# spike in params
total_ncells = 100
proportions = [0.3, 0.3, 0.4]

# bhc params
n_states = 8
alpha = 3
prob_cn_change = 0.9
bhc_incon = 2 # inconsistent score used for making clusters from bhc
bhc_depth = 2 

# naive clusering params
naive_method = "complete"
naive_metric = "cityblock"
naive_incon = 1.1
naive_depth = 2

# Params for testing threshold values
params = simulation.expand_grid({
  "transform":["log","none"], 
  "criterion": ["inconsistent"], "threshold": np.arange(0.025, 2, step=0.05)
})
params = pd.concat([params, 
  simulation.expand_grid({
    "transform":["log","none"], 
    "criterion": ["distance"], 
    "threshold": np.arange(3, 20, step=1)
  })
])

start = time.time()
bhc_linkage, bhc_root, bhc_cell_ids, matrix_data, measurement, variances = (
    cncluster.bayesian_cluster(cn_data, n_states=n_states, alpha=alpha, 
                               prob_cn_change=prob_cn_change,
                               clustering_id="copy", debug=True,
                               print_r=True)
)
print(f"{time.time()-start}s for BHC on {total_ncells} cells")

bhc_linkage, bhc_plot_data = simulation.get_plot_data(bhc_linkage)
lbhc_plot_data = bhc_plot_data.copy()
lbhc_plot_data[:,2] = np.log(bhc_plot_data[:,2]) 

bhc_clusters = sch.fcluster(lbhc_plot_data, 8, criterion="distance")
assert len(set(bhc_clusters)) > 1
cn_data = cncluster.prune_cluster(bhc_clusters, bhc_cell_ids, cn_data)
cn_data["origin_id_int"] = cn_data["origin_id"].factorize()[0]

fig = plt.figure(figsize=(10, 8))
bimatrix_data = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="bhc_cluster_id",
    linkage=lbhc_plot_data, origin_field_name="origin_id_int", raw=True, 
    flip=True)
clabels = get_mixture_labels(cn_data)
scores = skm.homogeneity_completeness_v_measure(
  clabels["origin_id_int"], clabels["bhc_cluster_id"]
)
print(f"homogeneity: {scores[0]}, completeness: {scores[1]}, "
      "v-measure: {scores[2]}")
