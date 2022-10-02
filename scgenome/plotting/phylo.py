import Bio.Phylo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import gridspec
import seaborn as sns

import scgenome.plotting.heatmap


def map_annotations_to_colors(annotation, cmap):
    num_colors = len(np.unique(annotation))

    color_map = {}
    for idx, value in enumerate(np.unique(annotation)):
        color_map[value] = cmap(float(idx) / float(max(1, num_colors - 1)))

    color_mat = []
    for value in annotation:
        color_mat.append(color_map[value])
        
    return color_mat, color_map


def plot_tree_cn(
        tree,
        adata,
        chrom_segments=True,
        layer_name=None,
        obs_annotation=None,
        obs_cmap=None,
        var_label=None,
        fig=None,
        raw=False,
        max_cn=13):
    """ Plot a tree aligned to a CN values matrix heatmap

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        phylogenetic tree
    adata : AnnData
        Copy number data, either genes or segments
    chrom_segments : bool, optional
        whether adata is genes or segments, by default True
    layer_name : str, optional
        layer to plot for copy number heatmap, by default None
    obs_annotation : str, optional
        column of adata.obs to annotate cells, by default None
    obs_cmap : matplotlib.colors.ListedColormap, optional
        color map for cell annotations, by default None
    var_label : str, optional
        column of adata.var to use for heatmap var labels, by default None use index
    fig : matplotlib.figure.Figure, optional
        existing figure to plot into, by default None
    raw : bool, optional
        raw plotting, no integer color map, by default False
    max_cn : int, optional
        clip cn at max value, by default 13
    """    
    
    if fig is None:
        fig = plt.figure(figsize=(16, 12), dpi=150)

    # Add phylogenetic ordering to anndata obs
    cell_ids = []
    for a in tree.get_terminals():
        cell_ids.append(a.name)

    adata.obs['phylo_order'] = None
    for idx, _ in adata.obs.iterrows():
        adata.obs.loc[idx, 'phylo_order'] = cell_ids.index(idx)

    gs = gridspec.GridSpec(1, 3, width_ratios=(0.4, 0.58, 0.02))

    # Plot phylogenetic tree
    ax = fig.add_subplot(gs[0, 0])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    Bio.Phylo.draw(tree, label_func=lambda a: '', axes=ax, do_show=False)

    # Plot copy number heatmap
    ax = fig.add_subplot(gs[0, 1])

    if chrom_segments:
        scgenome.plotting.heatmap.plot_cell_cn_matrix(
            adata,
            layer_name=layer_name,
            cell_order_fields=('phylo_order',),
            ax=ax,
            raw=raw,
            max_cn=max_cn,
        )

    else:
        if layer_name is not None:
            X = adata.layers[layer_name]
        else:
            X = adata.X

        # Order the cells according to the phylogeny
        cell_order_values = adata.obs[['phylo_order']].values.transpose()
        cell_ordering = np.lexsort(cell_order_values)

        X = X[cell_ordering, :]

        mat = adata[cell_ordering, :].to_df(layer=layer_name)
        if var_label is not None:
            mat = mat.rename(columns=adata.var[var_label])
        
        cbar_ax = fig.add_axes([0.45, 0.0, 0.2, 0.01])
        sns.heatmap(mat, ax=ax, cbar_ax=cbar_ax, cbar_kws={'orientation': 'horizontal'})

    ax.set_yticks([])

    # Plot obs annotation
    if obs_annotation is not None:
        ax = fig.add_subplot(gs[0, 2])

        color_mat, color_map = map_annotations_to_colors(adata.obs[obs_annotation], obs_cmap)
        color_mat = np.swapaxes(np.array([color_mat]), 0, 1)

        ax.imshow(color_mat[::-1, :, :], aspect='auto', origin='lower', interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])

        patches = []
        for label, color in color_map.items():
            patches.append(mpatches.Patch(color=color, label=label))
        ax.legend(handles=patches, loc='upper left', bbox_to_anchor=(1.05, 1.))

    plt.subplots_adjust(left=0.065, right=0.97, top=0.96, bottom=0.065, wspace=0.01)

    return fig
