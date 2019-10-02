# Script for selecting parameters for pruning tree
import pandas as pd
from os.path import join
import numpy as np
from scgenome import simulation, cncluster, cnplot
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib.pyplot as plt

arg = {
    'out_dir': "/Users/massoudmaher/data/mixture1",
    'linkage_fp': "/Users/massoudmaher/data/mixture1/bhc_plot_data.txt",
    'cn_data_fp': '/Users/massoudmaher/data/mixture1/cn_data.csv',
    'cell_ids_fp': '/Users/massoudmaher/data/mixture1/bhc_cell_ids.csv',
    'cnd_out_fn': 'bhc_clst_cn_data.csv',
    # What cutoff values/transformations to try
    'grid': {
        "transform": ["none"],
        "criterion": ["distance"],
        "threshold": np.arange(1e6, 3e6, step=1e5)
    },
    'selection': {
        "transform": "none",
        "criterion": "distance",
        "threshold": 1.2e6
    },
    'plot': True,
    'spike_in': True,
    'prefix': 'bhc'
}
cnd_out = ('sel_clst_cn_data.csv'if arg['cnd_out_fn'] is None
           else arg['cnd_out_fn'])
origin = 'origin_id_int' if arg['spike_in'] else None
cluster_field_name = arg['prefix'] + "_cluster_id"



linkage = np.loadtxt(arg['linkage_fp'])
llinkage = linkage.copy()
llinkage[:, 2] = np.log(llinkage[:, 2])

grid = simulation.expand_grid(arg['grid'])


def apply_fn(row):
    if row["transform"] == "log":
        df = llinkage
    else:
        df = linkage
    return sch.fcluster(df, row["threshold"], criterion=row["criterion"])
grid["fcluster"] = grid.apply(apply_fn, axis=1)
grid["num_clusters"] = grid["fcluster"].apply(lambda x: len(set(x)))

fig = plt.figure(figsize=(8, 8))
g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                  sharey=True, sharex=False).add_legend()
g = g.map(plt.scatter, "threshold", "num_clusters")
plt.savefig(join(arg['out_dir'], arg['prefix'] + "_cluster_thresh.png"))

if arg['selection'] is not None:
    if arg['cn_data_fp'] is None:
        raise ValueError('selection must be provided with cn_data_fp')
    sel = arg['selection']
    z = llinkage if sel['transform'] == 'log' else linkage
    sel_clst = sch.fcluster(z, sel['threshold'], criterion=sel['criterion'])
    if len(set(sel_clst)) <= 1:
        raise AssertionError("Only have 1 cluster")

    cn_data = pd.read_csv(arg['cn_data_fp'])
    cell_ids = pd.read_csv(arg['cell_ids_fp'], header=None)
    cell_ids = list(cell_ids.iloc[:, 1])
    print(cell_ids)
    print(sel_clst)
    print(len(sel_clst))
    cn_data = cncluster.prune_cluster(sel_clst, cell_ids, cn_data,
                                      cluster_field_name)

    cn_data.to_csv(join(arg['out_dir'], cnd_out))

if arg['selection'] is not None and arg['plot']:
    # NHC

    fig = plt.figure(figsize=(9.4, 9))
    bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, cn_data, "copy", cluster_field_name=cluster_field_name,
        linkage=z, origin_field_name=origin, raw=True,
        flip=True, cell_id_order=cell_ids)
    fig.savefig(join(arg['out_dir'], arg['prefix'] + "_re_heatmap.png"),
                bbox_inches='tight')

print(cn_data.shape)
