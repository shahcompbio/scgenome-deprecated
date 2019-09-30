# Script for selecting parameters for pruning tree
import pandas as pd
import os
import numpy as np
from scgenome import simulation, cncluster
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib.pyplot as plt

arg = {
    'out_dir': "/Users/massoudmaher/data/bhc_sc1935",
    'linkage_fp': "/Users/massoudmaher/data/bhc_sc1935/bhc_plot_data.txt",
    'llinkage_fp': "/Users/massoudmaher/data/bhc_sc1935/lbhc_plot_data.txt",
    'cn_data_fp': '/Users/massoudmaher/data/sc1935_clean_qc.csv',
    # What cutoff values/transformations to try
    'grid': {
        "transform": ["log", "none"],
        "criterion": ["distance"],
        "threshold": np.arange(3, 20, step=0.5)
    },
    'selection': {
        "transform": "log",
        "criterion": "distance",
        "threshold": 12.5
    }
}


linkage = np.loadtxt(arg['linkage_fp'])
llinkage = np.loadtxt(arg['llinkage_fp'])

grid = simulation.expand_grid(arg['grid'])


def apply_fn(row):
    if row["transform"] == "log":
        df = llinkage
    else:
        df = linkage
    return sch.fcluster(df, row["threshold"], criterion=row["criterion"])
grid["bhc_fcluster"] = grid.apply(apply_fn, axis=1)
grid["bhc_num_clusters"] = grid["bhc_fcluster"].apply(lambda x: len(set(x)))

fig = plt.figure(figsize=(8, 8))
g = sns.FacetGrid(data=grid, col="criterion", row="transform", size=6,
                  sharey=True, sharex=False).add_legend()
g = g.map(plt.scatter, "threshold", "bhc_num_clusters")
plt.savefig(os.path.join(arg['out_dir'], "cluster_thresh.pdf"))

if arg['selection'] is not None:
    if arg['cn_data_fp'] is None:
        raise ValueError('selection must be provided with cn_data_fp')
    sel = arg['selection']
    z = llinkage if sel['transform'] == 'log' else linkage
    sel_clst = sch.fcluster(z, sel['threshold'], criterion=sel['criterion'])
    if len(set(sel_clst)) <= 1:
        raise AssertionError("Only have 1 cluster")

    cn_data = pd.read_csv(arg['cn_data_fp'])
    cell_ids = pd.read_csv(os.path.join(arg['out_dir'], 'cell_ids.csv'),
                           header=None)
    cell_ids = list(cell_ids.iloc[:, 1])
    print(cell_ids)
    print(sel_clst)
    print(len(sel_clst))
    cn_data = cncluster.prune_cluster(sel_clst, cell_ids, cn_data)
    cn_data["origin_id_int"] = cn_data["origin_id"].factorize()[0]

    cn_data.to_csv(os.path.join(arg['out_dir'], 'clst_cn_data.csv'))



