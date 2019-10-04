# Do Naive clustering on a dataset
import pandas as pd
import numpy as np
from scgenome import cncluster, utils, simulation, cnplot
import scipy.cluster.hierarchy as sch
import os
import time
import matplotlib.pyplot as plt
import sklearn.metrics as skm

from scgenome.constants import VALUE_IDS

arg = {
    'out_dir': "/Users/massoudmaher/data/test_do_naive/",
    'cn_data_fp': "/Users/massoudmaher/data/"
                  "clean_sc_1935_1936_1937_cn_data_qc.csv",
    'seed': None,

    'n_cells': 200,  # Set to None if we want all cells or all bins
    'n_bin': None,

    # naive params
    'naive_method': 'complete',
    'naive_metric': 'euclidean',

    # spike in params (not always used)
    'sample_ids':  ['SC-1935', 'SC-1936', 'SC-1937'],
    'spike_in': True,
    'proportions': None, # Set to None for equal proportion of each sample
    'nhc_thresh': 70
}

pd.DataFrame(arg).to_csv(os.path.join(arg['out_dir'], "do_naive_arg.csv"))

if not os.path.exists(arg['out_dir']):
    print(f"{arg['out_dir']} does not exist, creating it")
    os.makedirs(arg['out_dir'])

print(f"Reading in {arg['cn_data_fp']}")
cn_data = pd.read_csv(arg['cn_data_fp'])

if arg['n_cells'] is not None:
    if arg['spike_in']:
        cn_data = utils.get_cn_data_submixture(cn_data, arg['n_cells'],
                                               arg['sample_ids'],
                                               proportions=arg['proportions'])
    else:
        cells = pd.Series(np.unique(cn_data["cell_id"]))
        keep_cells = cells.sample(arg['n_cells'], random_state=arg['seed'])
        cn_data = cn_data[cn_data["cell_id"].isin(keep_cells)]

if arg['n_bin'] is not None:
    end_val = np.unique(cn_data["end"])[arg['n_bin']]
    print(f"Reducing to {arg['n_bin']} bins, end: {end_val}")
    cn_data = cn_data[cn_data["end"] <= end_val]

if arg['n_cells'] is not None or arg['n_bin'] is not None:
    cn_data.to_csv(os.path.join(arg['out_dir'], "cn_data.csv"))

n_cell = np.unique(cn_data["cell_id"]).shape[0]
n_bin = np.unique(cn_data["end"]).shape[0]

print(f"cn_data.shape {cn_data.shape}")

matrix_data, measurement, cell_ids = (
    utils.cn_data_to_mat_data_ids(cn_data, data_id="copy",
                                  value_ids=VALUE_IDS))
cell_ids.to_csv(os.path.join(arg['out_dir'], 'cell_ids.csv'))

print(f"Doing NHC on {n_cell} cells, {n_bin} bins")
start = time.time()
naive_linkage = sch.linkage(np.nan_to_num(measurement),
                            method=arg['naive_method'],
                            metric=arg['naive_metric'])
print(f"{time.time()-start}s for NHC on {n_cell} cells, {n_bin} bins\n\n")
np.savetxt(os.path.join(arg['out_dir'], "naive_plot_data.txt"), naive_linkage)


print("Plotting")
# NHC
naive_clusters = sch.fcluster(naive_linkage, arg['nhc_thresh'],
                              criterion="distance")
assert len(set(naive_clusters)) > 1
cn_data = cncluster.prune_cluster(naive_clusters, cell_ids, cn_data,
                                  cluster_field_name="naive_cluster_id")

origin = 'origin_id_int' if arg['spike_in'] else None

fig = plt.figure(figsize=(10, 8))
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="naive_cluster_id",
    linkage=naive_linkage, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=cell_ids)
fig.savefig(os.path.join(arg['out_dir'], "naive_heatmap.png"),
            bbox_inches='tight')

