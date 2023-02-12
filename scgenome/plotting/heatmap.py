import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

from anndata import AnnData

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch
import matplotlib.cm

import scgenome.cnplot
import scgenome.refgenome
from . import cn_colors


def plot_cell_cn_matrix(
        adata: AnnData,
        layer_name='state',
        cell_order_fields=(),
        ax=None,
        raw=False,
        max_cn=13,
        vmin=None,
        vmax=None,
        show_cell_ids=False):
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
    vmin, vmax : float, optional
        for raw=True, vmin and vmax define the data range that the colormap covers, see `matplotlib.pyplot.imshow`
    show_cell_ids : bool, optional
        show cell ids on heatmap axis, by default False

    Returns
    -------
    dict
        Dictionary of plot and data elements

    Examples
    -------

    .. plot::
        :context: close-figs

        import scgenome
        adata = scgenome.datasets.OV2295_HMMCopy_reduced()
        scgenome.pl.plot_cell_cn_matrix(adata)

    """

    if ax is None:
        ax = plt.gca()

    # Order the chromosomes
    chr_start = adata.var.reset_index().merge(scgenome.refgenome.info.chromosome_info[['chr', 'chr_index']], how='left')
    if chr_start['chr_index'].isnull().any():
        chromosomes = adata.var['chr'].astype(str).values
        raise ValueError(f'mismatching chromosomes {chromosomes} and {scgenome.refgenome.info.chromosomes}')
    chr_start = chr_start[['start', 'chr_index']].values
    genome_ordering = np.lexsort(chr_start.transpose())

    # Order the cells if requested
    if len(cell_order_fields) > 0:
        cell_order_fields = reversed(list(cell_order_fields))
        cell_order_values = adata.obs[cell_order_fields].values.transpose()
        cell_ordering = np.lexsort(cell_order_values)

    else:
        cell_ordering = range(adata.shape[0])

    adata = adata[cell_ordering, genome_ordering]

    if layer_name is not None:
        X = adata.layers[layer_name].copy()
    else:
        X = adata.X.copy()

    X = np.nan_to_num(X, nan=0)

    if max_cn is not None:
        X[X > max_cn] = max_cn

    if not raw:
        X_colors = cn_colors.map_cn_colors(X)
        im = ax.imshow(X_colors, aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)

    else:
        cmap = matplotlib.cm.get_cmap('viridis')
        im = ax.imshow(X, aspect='auto', cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)

    mat_chrom_idxs = chr_start[genome_ordering][:, 1]
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [mat_chrom_idxs.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2
    ordered_mat_chrom_idxs = mat_chrom_idxs[np.where(np.array([1] + list(np.diff(mat_chrom_idxs))) != 0)]
    chrom_names = np.array(scgenome.refgenome.info.plot_chromosomes)[ordered_mat_chrom_idxs]

    ax.set_xticks(chrom_mids)
    ax.set_xticklabels(chrom_names, fontsize='6')

    if show_cell_ids:
        ax.set(yticks=range(len(adata.obs.index)))
        ax.set(yticklabels=adata.obs.index.values)
        ax.tick_params(axis='y', labelrotation=0)
    else:
        ax.set(yticks=[])
        ax.set(yticklabels=[])

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)

    return {
        'ax': ax,
        'im': im,
        'adata': adata,
    }


def map_catagorigal_colors(values, cmap_name=None):
    levels = np.unique(values)
    n_levels = len(levels)

    if cmap_name is None:
        if n_levels <= 10:
            cmap_name = 'tab10'
        elif n_levels <= 20:
            cmap_name = 'tab20'
        else:
            cmap_name = 'hsv'

    cmap = matplotlib.cm.get_cmap(cmap_name)

    level_colors = dict(zip(levels, cmap(np.linspace(0, 1, n_levels))))

    value_colors = np.zeros(values.shape + (4,))
    for l, c in level_colors.items():
        value_colors[values == l, :] = c

    return level_colors, value_colors


def _plot_categorical_annotation(values, ax, ax_legend, title):
    level_colors, value_colors = map_catagorigal_colors(values)

    im = ax.imshow(value_colors, aspect='auto', interpolation='none')

    levels = []
    patches = []
    for s, h in level_colors.items():
        levels.append(s)
        patches.append(Patch(facecolor=h, edgecolor=h))
    ncol = min(3, int(len(levels)**(1/2)))
    legend = ax_legend.legend(patches, levels, ncol=ncol,
        frameon=True, loc=2, bbox_to_anchor=(0., 1.),
        facecolor='white', edgecolor='white', fontsize='4',
        title=title, title_fontsize='6')

    ax.grid(False)
    ax.set_xticks([0.], [title], rotation=90, fontsize='6')
    ax.tick_params(axis='y', left=False, right=False)

    annotation_info = {}
    annotation_info['ax'] = ax
    annotation_info['im'] = im
    annotation_info['ax_legend'] = ax_legend
    annotation_info['legend'] = legend
    annotation_info['level_colors'] = level_colors
    annotation_info['value_colors'] = value_colors

    return annotation_info


def _plot_continuous_legend(ax_legend, im, title):
    ax_legend.grid(False)
    ax_legend.set_xticks([])
    ax_legend.set_yticks([])

    axins = ax_legend.inset_axes([0.5, 0.1, 0.05, 0.8])

    cbar = plt.colorbar(im, cax=axins)
    axins.set_title(title, fontsize='6')
    cbar.ax.tick_params(labelsize='4')

    annotation_info = {}
    annotation_info['ax_legend'] = ax_legend
    annotation_info['axins'] = axins
    annotation_info['cbar'] = cbar

    return annotation_info


def _plot_continuous_annotation(values, ax, ax_legend, title):

    im = ax.imshow(values, aspect='auto', interpolation='none', cmap='Reds')

    ax.grid(False)
    ax.set_xticks([0.], [title], rotation=90, fontsize='6')
    ax.tick_params(axis='y', left=False, right=False)

    annotation_info = _plot_continuous_legend(ax_legend, im, title)

    annotation_info['ax'] = ax
    annotation_info['im'] = im

    return annotation_info


def plot_cell_cn_matrix_fig(
        adata: AnnData,
        layer_name='state',
        cell_order_fields=None,
        annotation_fields=None,
        fig=None,
        raw=False,
        vmin=None,
        vmax=None,
        max_cn=13,
        show_cell_ids=False,
        show_subsets=False):
    """ Plot a copy number matrix

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    cell_order_fields : list, optional
        columns of obs on which to sort cells, by default None
    annotation_fields : list, optional
        column of obs to use as an annotation colorbar, by default 'cluster_id'
    fig : matplotlib.figure.Figure, optional
        existing figure to plot into, by default None
    raw : bool, optional
        raw plotting, no integer color map, by default False
    vmin, vmax : float, optional
        for raw=True, vmin and vmax define the data range that the colormap covers, see `matplotlib.pyplot.imshow`
    max_cn : int, optional
        clip cn at max value, by default 13
    show_cell_ids : bool, optional
        show cell ids on heatmap axis, by default False
    show_subsets : bool, optional
        show subset/superset categoricals to allow identification of cell sets

    Returns
    -------
    dict
        Dictionary of plot and data elements

    Examples
    -------

    .. plot::
        :context: close-figs

        import scgenome
        adata = scgenome.datasets.OV2295_HMMCopy_reduced()

        g = scgenome.pl.plot_cell_cn_matrix_fig(
            adata,
            cell_order_fields=['cell_order'],
            annotation_fields=['cluster_id', 'sample_id', 'quality'])

    """

    if fig is None:
        fig = plt.figure()

    if cell_order_fields is None:
        cell_order_fields = []

    if annotation_fields is None:
        annotation_fields = []

    # Account for additional annotation fields that will be added after plotting the matrix
    # when we are adding annotation fields to identify cells as per show_subsets
    num_annotations = len(annotation_fields)
    if show_subsets:
        num_annotations += 2

    fig_main, fig_legends = fig.subfigures(nrows=2, ncols=1, height_ratios=[5, 1], squeeze=True)
    fig_legends.patch.set_alpha(0.0)

    width_ratios = [1] + [0.005] + [0.02] * num_annotations

    axes = fig_main.subplots(
        nrows=1, ncols=len(width_ratios), sharey=True, width_ratios=width_ratios,
        squeeze=False, gridspec_kw=dict(hspace=0.02, wspace=0.02))[0]
    axes[1].set_axis_off()

    axes_legends = fig_legends.subplots(
        nrows=1, ncols=1+num_annotations, squeeze=False)[0]
    for ax in axes_legends:
        ax.set_axis_off()
        ax.set_alpha(0.0)
        ax.patch.set_alpha(0.0)

    ax = axes[0]
    ax_legend = axes_legends[0]
    g = plot_cell_cn_matrix(
        adata, layer_name=layer_name,
        cell_order_fields=cell_order_fields,
        ax=ax, raw=raw, vmin=vmin, vmax=vmax,
        max_cn=max_cn, show_cell_ids=show_cell_ids)

    adata = g['adata']
    im = g['im']

    if not raw:
        legend_info = {'ax_legend': ax_legend}
        legend_info['legend'] = cn_colors.cn_legend(ax_legend, title=layer_name)

    else:
        legend_info = _plot_continuous_legend(ax_legend, im, layer_name)

    if show_subsets:
        # Need to copy the adata to avoid modifying a view
        adata = adata.copy()
        adata.obs['subset'] = pd.Series(np.mod(np.floor_divide(range(adata.shape[0]), 40), 5), index=adata.obs.index, dtype='category')
        adata.obs['superset'] = pd.Series(np.floor_divide(range(adata.shape[0]), 200), index=adata.obs.index, dtype='category')
        annotation_fields.append('superset')
        annotation_fields.append('subset')

    annotation_info = {}

    for ax, ax_legend, annotation_field in zip(axes[2:], axes_legends[1:], annotation_fields):
        if adata.obs[annotation_field].dtype.name in ('category', 'object'):
            values = adata.obs[[annotation_field]].values
            annotation_info[annotation_field] = _plot_categorical_annotation(values, ax, ax_legend, annotation_field)

        else:
            values = adata.obs[[annotation_field]].values
            annotation_info[annotation_field] = _plot_continuous_annotation(values, ax, ax_legend, annotation_field)

    return {
        'fig': fig,
        'axes': axes,
        'adata': adata,
        'im': im,
        'legend_info': legend_info,
        'annotation_info': annotation_info,
    }
