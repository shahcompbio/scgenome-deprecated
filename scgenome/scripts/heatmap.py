import matplotlib.pyplot as plt
import os
from scgenome import cnplot
import pandas as pd
import numpy as np

arg = {
    'out_dir': "/Users/massoudmaher/data/bhc_sc1935",
    'cn_data_fp': '/Users/massoudmaher/data/sc1935_clean_qc.csv',
    # What cutoff values/transformations to try
    'selection': {
        "transform": "log",
        "criterion": "distance",
        "threshold": 12.5
    },
    'out_name': 'log_dist_12p5.pdf',
    'spike_in': False
}

cn_data = pd.read_csv(os.path.join(arg['out_dir'], 'clst_cn_data.csv'))
linkage = np.loadtxt(os.path.join(arg['out_dir'], 'lbhc_plot_data.txt'))
cell_ids = pd.read_csv(os.path.join(arg['out_dir'], 'cell_ids.csv'),
                       header=None)
cell_ids = cell_ids.iloc[:, 1]

fig = plt.figure(figsize=(10, 8))
origin = "origin_id_int" if arg['spike_in'] else None
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name="bhc_cluster_id",
    linkage=linkage, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=cell_ids)
fig.savefig(os.path.join(arg['out_dir'], arg['out_name']), bbox_inches='tight')


