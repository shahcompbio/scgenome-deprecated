import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

from anndata import AnnData

import scgenome.cnplot
import scgenome.refgenome


def plot_cell_cn_matrix(adata: AnnData, layer_name='state', cell_order_fields=None, ax=None, raw=False, max_cn=13):
    """ Plot a copy number matrix

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    cell_order_fields : list, optional
        columns of obs on which to sort cells, by default None
    ax : matplotlib.axes.Axes, optional
        existing axis to plot into, by default None
    raw : bool, optional
        raw plotting, no integer color map, by default False
    max_cn : int, optional
        clip cn at max value, by default 13

    TODO: missing return
    """    

    if ax is None:
        ax = plt.gca()

    if layer_name is not None:
        X = adata.layers[layer_name].copy()
    else:
        X = adata.X.copy()

    X = np.nan_to_num(X, nan=0)

    if max_cn is not None:
        X[X > max_cn] = max_cn

    # Order the chromosomes
    chr_start = adata.var.reset_index().merge(scgenome.refgenome.info.chrom_idxs, how='left')[['start', 'chr_index']].values
    genome_ordering = np.lexsort(chr_start.transpose())

    # Order the cells if requested
    if cell_order_fields is not None:
        cell_order_fields = reversed(list(cell_order_fields))
        cell_order_values = adata.obs[cell_order_fields].values.transpose()
        cell_ordering = np.lexsort(cell_order_values)

    else:
        cell_ordering = range(X.shape[0])

    X = X[cell_ordering, :][:, genome_ordering]

    cmap = None
    if not raw:
        cmap = scgenome.cnplot.get_cn_cmap(X)

    im = ax.imshow(X, aspect='auto', cmap=cmap, interpolation='none')

    mat_chrom_idxs = chr_start[genome_ordering][:, 1]
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [mat_chrom_idxs.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2
    ordered_mat_chrom_idxs = mat_chrom_idxs[np.where(np.array([1] + list(np.diff(mat_chrom_idxs))) != 0)]
    chrom_names = np.array(scgenome.refgenome.info.chromosomes)[ordered_mat_chrom_idxs]

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=chrom_names)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return ax


def plot_cell_cn_matrix_clusters_fig(
        adata: AnnData,
        layer_name='state',
        cell_order_fields=None,
        annotation_field='cluster_id',
        fig=None,
        raw=False,
        max_cn=13):
    """ Plot a copy number matrix

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    cell_order_fields : list, optional
        columns of obs on which to sort cells, by default None
    annotation_field : str
        column of obs to use as an annotation colorbar, by default 'cluster_id'
    fig : matplotlib.figure.Figure, optional
        existing figure to plot into, by default None
    raw : bool, optional
        raw plotting, no integer color map, by default False
    max_cn : int, optional
        clip cn at max value, by default 13
    """    

    if fig is None:
        fig = plt.figure()

    ax = fig.add_axes([0.1,0.0,0.9,1.])
    plot_cell_cn_matrix(
        adata, layer_name=layer_name,
        cell_order_fields=cell_order_fields,
        ax=ax, raw=raw, max_cn=max_cn)

    cluster_ids = adata.obs.sort_values(cell_order_fields)[annotation_field].values
    color_mat = scgenome.cncluster.get_cluster_colors(cluster_ids)

    ax = fig.add_axes([0.0,0.0,0.05,1.])
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower', interpolation='none')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    return fig
