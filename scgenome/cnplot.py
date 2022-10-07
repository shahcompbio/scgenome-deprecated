import logging
import matplotlib
import seaborn
import random
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dst
from sklearn.decomposition import PCA
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


def _secondary_clustering(data):
    D = dst.squareform(dst.pdist(data.T, 'cityblock'))
    Y = sch.linkage(D, method='complete')
    Z = sch.dendrogram(Y, color_threshold=-1, no_plot=True)
    idx = np.array(Z['leaves'])
    ordering = np.zeros(idx.shape[0], dtype=int)
    ordering[idx] = np.arange(idx.shape[0])
    return ordering


def plot_clustered_cell_cn_matrix(ax, cn_data, cn_field_name, cluster_field_name='cluster_id', secondary_field_name=None, raw=False, max_cn=13, cmap=None):
    plot_data = cn_data.merge(refgenome.info.chromosome_info[['chr', 'chr_index']])

    if secondary_field_name is not None:
        plot_data = plot_data.set_index(['chr_index', 'start', 'cell_id', cluster_field_name])
        plot_data = plot_data[[secondary_field_name, cn_field_name]]
        plot_data = plot_data.unstack(level=['cell_id', cluster_field_name]).fillna(0)
        ordering = plot_data[secondary_field_name].values[0]

        plot_data = cn_data.merge(refgenome.info.chromosome_info[['chr', 'chr_index']])
        plot_data = plot_data.set_index(['chr_index', 'start', 'cell_id', cluster_field_name])[cn_field_name].unstack(level=['cell_id', cluster_field_name]).fillna(0)
    else:
        plot_data = plot_data.set_index(['chr_index', 'start', 'cell_id', cluster_field_name])[cn_field_name].unstack(level=['cell_id', cluster_field_name]).fillna(0)
        ordering = _secondary_clustering(plot_data.values)
    ordering = pd.Series(ordering, index=plot_data.columns, name='cell_order')
    plot_data = plot_data.T.set_index(ordering, append=True).T

    plot_data = plot_data.sort_index(axis=1, level=[1, 2])
    if max_cn is not None:
        plot_data[plot_data > max_cn] = max_cn

    mat_chrom_idxs = plot_data.index.get_level_values(0).values
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [plot_data.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2
    ordered_mat_chrom_idxs = mat_chrom_idxs[np.where(np.array([1] + list(np.diff(mat_chrom_idxs))) != 0)]
    chrom_names = np.array(refgenome.info.plot_chromosomes)[ordered_mat_chrom_idxs]

    mat_cluster_ids = plot_data.columns.get_level_values(1).values
    cluster_boundaries = np.array([0] + list(np.where(mat_cluster_ids[1:] != mat_cluster_ids[:-1])[0]) + [plot_data.shape[1] - 1])
    cluster_sizes = cluster_boundaries[1:] - cluster_boundaries[:-1]
    cluster_mids = cluster_boundaries[:-1] + cluster_sizes / 2

    if not raw and cmap is None:
        cmap = get_cn_cmap(plot_data.values)

    im = ax.imshow(plot_data.astype(float).T, aspect='auto', cmap=cmap, interpolation='none')

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=chrom_names)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return plot_data


def plot_clustered_cell_cn_matrix_figure(fig, cn_data, cn_field_name, cluster_field_name='cluster_id', secondary_field_name=None, raw=False, max_cn=13, cmap=None):
    ax = fig.add_axes([0.1,0.0,0.9,1.])
    plot_data = plot_clustered_cell_cn_matrix(ax, cn_data, cn_field_name, cluster_field_name=cluster_field_name, secondary_field_name=secondary_field_name,raw=raw, max_cn=max_cn, cmap=cmap)

    cluster_ids = plot_data.columns.get_level_values(1).values
    color_mat = cncluster.get_cluster_colors(cluster_ids)

    ax = fig.add_axes([0.0,0.0,0.05,1.])
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower', interpolation='none')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    return plot_data


def plot_cell_cn_profile(ax, cn_data, value_field_name, cn_field_name=None, max_cn=13, chromosome=None, s=5, squashy=False, rawy=False, cmap=None):
    """ Plot copy number profile on a genome axis

    Args:
        ax: matplotlib axis
        cn_data: copy number table
        value_field_name: column in cn_data to use for the y axis value
    
    Kwargs:
        cn_field_name: state column to color scatter points
        max_cn: max copy number for y axis
        chromosome: single chromosome plot
        s: size of scatter points
        squashy: compress y axis
        rawy: raw data on y axis

    The cn_data table should have the following columns (in addition to value_field_name and
    optionally cn_field_name):
        - chr
        - start
        - end
    """
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start', 'chromosome_end']].copy()
    chromosome_info['chr'] = pd.Categorical(chromosome_info['chr'], categories=cn_data['chr'].cat.categories)
    plot_data = cn_data.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    squash_coeff = 0.15
    squash_f = lambda a: np.tanh(squash_coeff * a)
    if squashy:
        plot_data[value_field_name] = squash_f(plot_data[value_field_name])

    if cn_field_name is not None:
        if cmap is not None:
            ax.scatter(
                plot_data['start'], plot_data[value_field_name],
                c=plot_data[cn_field_name], s=s,
                cmap=cmap,
            )
        else:
            ax.scatter(
                plot_data['start'], plot_data[value_field_name],
                c=plot_data[cn_field_name], s=s,
                cmap=get_cn_cmap(plot_data[cn_field_name].astype(int).values),
            )
    else:
        ax.scatter(
            plot_data['start'], plot_data[value_field_name], s=s,
        )

    if chromosome is not None:
        chromosome_length = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_length']
        chromosome_start = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_start']
        chromosome_end = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_end']
        xticks = np.arange(0, chromosome_length, 2e7)
        xticklabels = ['{0:d}M'.format(int(x / 1e6)) for x in xticks]
        xminorticks = np.arange(0, chromosome_length, 1e6)
        ax.set_xlabel(f'chromosome {chromosome}')
        ax.set_xticks(xticks + chromosome_start)
        ax.set_xticklabels(xticklabels)
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xminorticks + chromosome_start))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim((chromosome_start, chromosome_end))

    else:
        ax.set_xlim((-0.5, refgenome.info.chromosome_info['chromosome_end'].max()))
        ax.set_xlabel('chromosome')
        ax.set_xticks([0] + list(refgenome.info.chromosome_info['chromosome_end'].values))
        ax.set_xticklabels([])
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_info['chromosome_mid']))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.plot_chromosomes))

    if squashy and not rawy:
        yticks = np.array([0, 2, 4, 7, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        ax.set_yticks(yticks_squashed)
        ax.set_yticklabels(ytick_labels)
        ax.set_ylim((-0.01, 1.01))
        ax.spines['left'].set_bounds(0, 1)
    elif not rawy:
        ax.set_ylim((-0.05*max_cn, max_cn))
        ax.set_yticks(range(0, int(max_cn) + 1))
        ax.spines['left'].set_bounds(0, max_cn)

    if chromosome is not None:
        seaborn.despine(ax=ax, offset=10, trim=False)
    else:
        seaborn.despine(ax=ax, offset=10, trim=True)

    return chromosome_info


def plot_breakends(ax, breakends, lw=0.5):
    """ Plot breakpoint flags and arcs on a genome axis

    Args:
        ax: matplotlib axis
        breakends: breakend table
    
    Kwargs:
        lw: line width for arcs

    The breakends table should have the following columns:
        - prediction_id
        - chromosome
        - position
        - strand

    """
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start']].copy()
    chromosome_info['chromosome'] = chromosome_info['chr']

    plot_data = breakends.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['plot_position'] = plot_data['position'] + plot_data['chromosome_start']

    xlim = ax.get_xlim()
    xrng = xlim[1] - xlim[0]

    ylim = ax.get_ylim()
    yrng = ylim[1] - ylim[0]

    arrow_length = 0.01 * xrng

    head_width = yrng / 30.
    head_length = xrng / 100.

    prediction_ids = list(plot_data['prediction_id'].unique())
    random.shuffle(prediction_ids)
    color_palette = seaborn.color_palette('hls', len(prediction_ids))
    prediction_colors = dict(zip(prediction_ids, color_palette))

    for idx in plot_data.index:
        prediction_id = plot_data.loc[idx, 'prediction_id']
        x = plot_data.loc[idx, 'plot_position']
        strand = plot_data.loc[idx, 'strand']
        y = np.random.uniform(ylim[0] + 0.75*yrng, ylim[0] + 0.95*yrng)
        offset = (arrow_length, -arrow_length)[strand == '+']
        color = prediction_colors[prediction_id]
        ax.arrow(x, y, offset, 0, color=color, lw=lw, alpha=1.0, head_width=head_width, head_length=head_length)
        ax.plot([x, x], [ylim[0], ylim[1]], color=color, lw=lw, ls='-', alpha=1.0)

    for prediction_id, pred_breakends in plot_data.groupby('prediction_id'):
        if len(pred_breakends.index) != 2:
            continue
        pos1 = pred_breakends['plot_position'].min()
        pos2 = pred_breakends['plot_position'].max()
        posmid = (pos1 + pos2) / 2.
        height = 0.5 * yrng * (pos2 - pos1) / xrng
        ypos = ylim[1]
        yposmid = ypos + height
        spl = scipy.interpolate.make_interp_spline([pos1, posmid, pos2], [ypos, yposmid, ypos], k=2)
        pos = np.linspace(pos1, pos2, 100)
        color = prediction_colors[prediction_id]
        ax.plot(pos, spl(pos), ls='-', lw=lw, color=color)

    ax.set_ylim((ylim[0], ylim[1] + 5))


def plot_cluster_cn_matrix(fig, cn_data, cn_field_name, cluster_field_name='cluster_id'):
    plot_data = cn_data.merge(refgenome.info.chromosome_info[['chr', 'chr_index']])
    plot_data = plot_data.groupby(['chr_index', 'start', cluster_field_name], observed=True)[cn_field_name].median().astype(int)
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
    chromosomes = refgenome.info.chromosome_info.set_index('chr_index').loc[np.unique(mat_chrom_idxs), 'chr_plot']

    ax = fig.add_axes([0.125,1.,0.875,1.])
    im = ax.imshow(plot_data.T, aspect='auto', cmap=get_cn_cmap(plot_data.values), interpolation='none')

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=chromosomes)
    ax.set(yticks=list(range(len(plot_data.columns.values))))
    ax.set(yticklabels=plot_data.columns.values)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return plot_data.columns.values


def compute_pca_loadings(cn_data):
    """ Compute the first n components of a PCA
    """
    cn_matrix = (
        cn_data
        .set_index(['cell_id', 'chr', 'start', 'end'])['copy']
        .unstack(level=[1, 2, 3]))

    num_null = cn_matrix.isnull().sum(axis=1)
    cn_matrix = cn_matrix[num_null <= 800]
    cn_matrix = cn_matrix.dropna(axis='columns')

    pca = PCA()
    pca.fit(cn_matrix.values)

    components = pd.DataFrame(
        pca.components_,
        columns=cn_matrix.columns)

    return components


def plot_pca_components(cn_data, n_components=4, plots_prefix=None):
    """ Plot the first n components of a PCA
    """
    components = compute_pca_loadings(cn_data)

    fig = plt.figure(figsize=(20, 4 * n_components))
    for idx in range(4):
        ax = fig.add_subplot(n_components, 1, idx + 1)
        plot_data = components.iloc[idx].T.rename('component').reset_index()
        plot_cell_cn_profile(
            ax, plot_data, 'component', rawy=True)
        ax.set_ylabel(f'PCA {idx+1}')

    if plots_prefix is not None:
        fig.savefig(plots_prefix + 'pca_components.pdf', bbox_inches='tight')

