import logging
import matplotlib
import seaborn
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dst
from matplotlib.colors import ListedColormap

from scgenome import refgenome
from scgenome import utils
from scgenome import cncluster


def hex_to_rgb(h):
    if h is None:
        return np.array((0, 0, 0))
    h = h.lstrip('#')
    return np.array(tuple(np.uint8(int(h[i:i+2], 16)) for i in (0, 2 ,4)))


color_reference = {0:'#3182BD', 1:'#9ECAE1', 2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33', 6:'#B30000', 7:'#980043', 8:'#DD1C77', 9:'#DF65B0', 10:'#C994C7', 11:'#D4B9DA'}


def get_cn_cmap(cn_data):
    min_cn = int(cn_data.min())
    max_cn = int(cn_data.max())
    assert min_cn - cn_data.min() == 0
    assert max_cn - cn_data.max() == 0
    color_list = []
    for cn in range(min_cn, max_cn+1):
        if cn > max(color_reference.keys()):
            cn = max(color_reference.keys())
        color_list.append(color_reference[cn])
    return ListedColormap(color_list)


def plot_cbar(ax):
    ax.imshow(np.array([np.arange(len(color_reference))]).T[::-1], cmap=cmap, aspect=1)
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(color_reference)))
    ax.set_yticklabels(np.arange(len(color_reference))[::-1])


def _secondary_clustering(data, linkage=None):
    if linkage is None:
        D = dst.squareform(dst.pdist(data.T, 'cityblock'))
        Y = sch.linkage(D, method='complete')
        Z = sch.dendrogram(Y, color_threshold=-1, no_plot=True)
        idx = np.array(Z['leaves'])
    else:
        Z = sch.dendrogram(linkage, color_threshold=-1, no_plot=True)
        idx = np.array(linkage)
    ordering = np.zeros(idx.shape[0], dtype=int)
    ordering[idx] = np.arange(idx.shape[0])
    return ordering


def plot_clustered_cell_cn_matrix(ax, cn_data, cn_field_name,
                                  cluster_field_name='cluster_id',
                                  raw=False, max_cn=13, linkage=None):
    plot_data = cn_data.merge(utils.chrom_idxs)
    plot_data = plot_data.set_index(['chr_index', 'start', 'cell_id', cluster_field_name])[cn_field_name].unstack(level=['cell_id', cluster_field_name]).fillna(0)

    ordering = _secondary_clustering(plot_data.values, linkage)
    ordering = pd.Series(ordering, index=plot_data.columns, name='cell_order')
    plot_data = plot_data.T.set_index(ordering, append=True).T

    plot_data = plot_data.sort_index(axis=1, level=[1, 2])
    if max_cn is not None:
        plot_data[plot_data > max_cn] = max_cn

    mat_chrom_idxs = plot_data.index.get_level_values(0).values
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [plot_data.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2

    mat_cluster_ids = plot_data.columns.get_level_values(1).values
    cluster_boundaries = np.array([0] + list(np.where(mat_cluster_ids[1:] != mat_cluster_ids[:-1])[0]) + [plot_data.shape[1] - 1])
    cluster_sizes = cluster_boundaries[1:] - cluster_boundaries[:-1]
    cluster_mids = cluster_boundaries[:-1] + cluster_sizes / 2

    cmap = None
    if not raw:
        cmap = get_cn_cmap(plot_data.values)

    im = ax.imshow(plot_data.T, aspect='auto', cmap=cmap)

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=utils.chrom_names)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return plot_data


def plot_clustered_cell_cn_matrix_figure(fig, cn_data, cn_field_name,
                                         cluster_field_name='cluster_id',
                                         raw=False, max_cn=13,
                                         linkage=None, origin_field_name=None):
    ax = fig.add_axes([0.1,0.0,0.9,1.])
    plot_data = plot_clustered_cell_cn_matrix(
        ax, cn_data, cn_field_name, cluster_field_name=cluster_field_name,
        raw=raw, max_cn=max_cn, linkage=linkage)

    cluster_ids = plot_data.columns.get_level_values(1).values
    color_mat = cncluster.get_cluster_colors(cluster_ids)

    if linkage is None:
        ax = fig.add_axes([0.0, 0.0, 0.05, 1.])
    else:
        ax = fig.add_axes([1.03, 0, 0.05, 1.])
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    return plot_data


def plot_cell_cn_profile(ax, cn_data, value_field_name, cn_field_name, max_cn=13, chromosome=None, s=5):
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start', 'chromosome_end']].copy()
    chromosome_info['chr'] = pd.Categorical(chromosome_info['chr'], categories=cn_data['chr'].cat.categories)
    plot_data = cn_data.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    ax.scatter(
        plot_data['start'], plot_data[value_field_name],
        c=plot_data[cn_field_name], s=s,
        cmap=get_cn_cmap(plot_data[cn_field_name].values),
    )

    ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.set_ylim((-0.5, max_cn))
    ax.set_yticks(np.arange(0, max_cn, 2))
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.chromosomes))

    seaborn.despine(ax=ax, offset=10, trim=True)

    return chromosome_info


def plot_cluster_cn_matrix(fig, cn_data, cn_field_name, cluster_field_name='cluster_id'):
    plot_data = cn_data.merge(utils.chrom_idxs)
    plot_data = plot_data.groupby(['chr_index', 'start', cluster_field_name])[cn_field_name].median().astype(int)
    plot_data = plot_data.unstack(level=[cluster_field_name]).fillna(0)
    plot_data = plot_data.sort_index(axis=1, level=1)

    logging.info(f'matrix with size {plot_data.T.shape}')

    ax = fig.add_axes([0.0,1.,0.1,1.])
    if plot_data.shape[1] > 2:
        D = dst.squareform(dst.pdist(plot_data.T, 'cityblock'))
        Y = sch.linkage(D, method='complete')
        Z = sch.dendrogram(Y, color_threshold=-1, orientation='left')
        idx = Z['leaves'][::-1]
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plot_data = plot_data.iloc[:, idx]

    mat_chrom_idxs = plot_data.index.get_level_values(0).values
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [plot_data.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2

    ax = fig.add_axes([0.125,1.,0.875,1.])
    im = ax.imshow(plot_data.T, aspect='auto', cmap=get_cn_cmap(plot_data.values))

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=utils.chrom_names)
    ax.set(yticks=list(range(len(plot_data.columns.values))))
    ax.set(yticklabels=plot_data.columns.values)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return plot_data.columns.values


def plot_bhc(measurement, linkage, cell_labels):
    return sns.clustermap(measurement, col_cluster=False, row_linkage=linkage,
                          yticklabels=cell_labels, cmap="YlGnBu")
