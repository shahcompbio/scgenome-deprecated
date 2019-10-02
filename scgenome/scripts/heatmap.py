import matplotlib.pyplot as plt
from os.path import join
from scgenome import cnplot
import pandas as pd
import numpy as np

arg = {
    'out_dir': "/Users/massoudmaher/data/mixture1",
    'cn_data_fp': '/Users/massoudmaher/data/mixture1/cn_data.csv',
    'out_name': 'bhc_new.png',
    'spike_in': True,
    'hierarchical': True,
    'cluster_field_name': 'bhc_cluster_id',
    'dummy_linkage': False
}

cn_data = pd.read_csv(join(arg['out_dir'], 'bhc_clst_cn_data.csv'))
cell_ids = pd.read_csv(join(arg['out_dir'], 'bhc_cell_ids.csv'), header=None)

##
cn_data = cn_data.iloc[:, 2:]
cell_ids = cell_ids.iloc[:, 1]
##

if arg['hierarchical']:
    linkage = np.loadtxt(join(arg['out_dir'], 'lbhc_plot_data.txt'))
else:
    linkage = None
    clusters = pd.read_csv(join(arg['out_dir'], 'umap_hdb_clust.csv'))
    cn_data = cn_data.merge(clusters)

print(cell_ids.head())

fig = plt.figure(figsize=(9.4, 9))
origin = "origin_id_int" if arg['spike_in'] else None
bimatrix_data, ps = cnplot.plot_clustered_cell_cn_matrix_figure(
    fig, cn_data, "copy", cluster_field_name=arg['cluster_field_name'],
    linkage=linkage, origin_field_name=origin, raw=True,
    flip=True, cell_id_order=cell_ids, dummy_linkage=arg['dummy_linkage'])
fig.savefig(join(arg['out_dir'], arg['out_name']), bbox_inches='tight')

