import Bio.Phylo
import colors
import numpy as np
import pandas as pd
import refgenome
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

"""
data is dataframe

format:
                         CELL1             CELL2
1:1-500000             1                  2
1:500000-1000000      3                  4



cell_annotation_data

cell_id label1  label2
CELL1   A   C
CELL2   B   C

bin_annotation_data

bin                label1  label2
1:1-500000         1       2
1:500000-1000000   3       4


usage:
df = pd.read_csv("small_dataset.csv.gz")
df["bin"] = df.apply(lambda row: f"{row['chr']}:{row['start']}-{row['end']}", axis=1)
df = df.pivot(
    index='cell_id',
    columns='bin',
    values='state')
hmap = HeatMap(df)
plt.savefig('out.png')

df = pd.read_csv("small_dataset.csv.gz")
cellids = list(df['cell_id'].unique())
annotations = pd.DataFrame(
    [(v, True if i < 3 else False, False if i < 2 else True, i) for i, v in enumerate(cellids)],
    columns=['cell_id', 'contaminated', 'control', 'num_reads']
)
annotations = annotations[['cell_id', 'contaminated', 'control']]

df["bin"] = df.apply(lambda row: f"{row['chr']}:{row['start']}-{row['end']}", axis=1)
df = df.pivot(
    index='cell_id',
    columns='bin',
    values='state')
bin_annotations = [
    [bin, (1 if i < 3000 else 2), (True if i < 3000 else False), ('gneg' if i < 3000 else 'gpos')] for i, bin in
    enumerate(df.columns.values)
]
bin_annotations = pd.DataFrame(bin_annotations, columns=['bin', 'val1', 'val2', 'cyto_band_giemsa_stain'])

hmap = AnnotatedHeatMap(df, cell_annotation_data=annotations, bin_annotation_data=bin_annotations,
                        cell_order=cellids[::-1], subset=True)
hmap.view_by_cell_annotation('subset', 1)
hmap.show_cell_labels()
hmap.hide_cell_labels()
hmap.reset_cells()
plt.savefig('out.png')



"""


class HeatMap(object):
    def __init__(
            self,
            data,
            ax=None,
            cmap=None,
            cell_order=None,
    ):

        self.data, self.cell_order, self.bins = self._reformat_data(data, cell_order)

        self.cmap = cmap

        self.ax = plt.gca() if ax is None else ax

        self.plot_heatmap()

    def _reformat_data(self, data, cell_order):
        data = data.fillna(0)

        if cell_order is None:
            cell_order = list(data.index)

        data.index = pd.Categorical(data.index, cell_order)
        data = data.sort_index()

        bins = self.get_sorted_bins(data)
        data = data[bins]

        return data, cell_order, bins

    def get_sorted_bins(self, data):

        bins = data.columns

        bins = bins.drop_duplicates()

        bins = bins.str.split(':|-', expand=True).to_frame()

        bins.columns = ['chr', 'start', 'end']

        bins['start'] = bins['start'].astype(int)
        bins['end'] = bins['end'].astype(int)

        # just the sort order in here, dont need to change for mouse.
        # might need to extend if reference has more genomes than human
        chromosomes = refgenome.chromosomes()
        if bins.iloc[0][0].startswith('chr'):
            chromosomes = ['chr' + v for v in chromosomes]
        bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

        bins = bins.sort_values(['chr', 'start', 'end'])

        bins = bins.apply(lambda row: f"{row['chr']}:{row['start']}-{row['end']}", axis=1)

        return bins

    def show_cell_labels(self):
        self.ax.set(yticks=range(len(self.cell_order)))
        self.ax.set(yticklabels=self.cell_order)
        self.ax.tick_params(axis='y', labelrotation=0)
        return self.ax

    def hide_cell_labels(self):
        self.ax.set(yticks=[])
        self.ax.set(yticklabels=[])
        return self.ax

    def plot_heatmap(self):

        mat = np.array(self.data)
        mat = np.nan_to_num(mat, nan=0)

        X_colors = colors.map_cn_colors_to_matrix(mat)

        self.ax.imshow(X_colors, aspect='auto', interpolation='none')

        chridxs = np.array(self.bins.str.split(':').str[0])
        mat_chrom_idxs = np.array(sorted(np.unique(chridxs, return_inverse=True)[1]))

        chrom_boundaries = np.array(
            [0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [mat_chrom_idxs.shape[0] - 1]
        )

        chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
        chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2
        ordered_mat_chrom_idxs = mat_chrom_idxs[np.where(np.array([1] + list(np.diff(mat_chrom_idxs))) != 0)]
        chrom_names = np.array(refgenome.plot_chromosomes())[ordered_mat_chrom_idxs]

        self.ax.set_xticks(chrom_mids)
        self.ax.set_xticklabels(chrom_names, fontsize='6')
        self.ax.set_xlabel('chromosome', fontsize=8)
        for val in chrom_boundaries[:-1]:
            self.ax.axvline(x=val, linewidth=0.5, color='black', zorder=100)

        self.hide_cell_labels()


class AnnotatedHeatMap(object):
    def __init__(
            self,
            data,
            cell_annotation_data=None,
            bin_annotation_data=None,
            cmap=None,
            cell_order=None,
            phylogenetic_tree=None,
            fig=None,
            subset=False
    ):

        self.subset = subset
        self.cmap = cmap
        self.tree = phylogenetic_tree

        if phylogenetic_tree is not None and cell_order is not None:
            raise Exception('both tree and cell_order provided. please provide only one')

        if cell_annotation_data is not None:
            self.cell_annotation_data = cell_annotation_data.set_index('cell_id')
            self.num_cell_annotations = len(self.cell_annotation_data.columns)
        else:
            self.cell_annotation_data = None
            self.num_cell_annotations = 0

        if bin_annotation_data is not None:
            self.bin_annotation_data = bin_annotation_data.set_index('bin')
            self.num_bin_annotations = len(self.bin_annotation_data.columns)
        else:
            self.bin_annotation_data = None
            self.num_bin_annotations = 0

        self.fig = plt.figure(figsize=(10, 10)) if fig is None else fig
        main_fig, legend_fig = self.setup_primary_fig(self.fig)
        self.tree_ax, self.heatmap_ax, self.cell_ann_axes, self.bin_ann_axes, self.subset_ax = self.setup_heatmap_axes(
            main_fig, disable_axis=True)
        self.heatmap_ax_legend, self.cell_ann_axes_legends, self.bin_ann_axes_legends, self.subset_ax_legend = self.setup_legend_axes(
            legend_fig, disable_axis=True)

        self.heatmap_ax.set_axis_on()

        for cell_annotation in self.cell_ann_axes:
            cell_annotation.set_axis_on()

        for bin_annotation in self.bin_ann_axes:
            bin_annotation.set_axis_on()

        self.heatmap = HeatMap(data, ax=self.heatmap_ax, cell_order=cell_order, cmap=cmap)

        self.cell_order = self.heatmap.cell_order
        self.bins = self.heatmap.bins

        self.plot_cell_annotations()

        self.plot_bin_annotations()

        self.plot_phylogenetic_tree()

        if self.subset:
            self.subset_ax.set_axis_on()
            self.add_subset_annotation()

        self.full_ylims = self.heatmap_ax.get_ylim()

    def setup_heatmap_axes(self, fig, disable_axis=False):

        num_cell_annotations = self.num_cell_annotations + 1 if self.subset else self.num_cell_annotations

        width_ratios = [1.0] + [0.02] * num_cell_annotations
        if self.tree:
            width_ratios = [0.5] + width_ratios

        height_ratios = [0.02] * self.num_bin_annotations + [1]

        axes_main = fig.subplots(
            nrows=len(height_ratios), ncols=len(width_ratios),
            width_ratios=width_ratios, height_ratios=height_ratios,
            squeeze=False, gridspec_kw=dict(hspace=0.02, wspace=0.02)
        )
        if disable_axis:
            for ax_row in axes_main:
                for ax in ax_row:
                    ax.set_axis_off()
                    ax.set_yticks([])
                    ax.set_xticks([])
        # first few rows are annotations, last one is heatmap and tree
        tree_ax = axes_main[-1][0] if self.tree is not None else None
        heatmap_ax = axes_main[-1][0] if self.tree is None else axes_main[-1][1]
        bin_ann_axes = [v[0] for v in axes_main] if self.tree is None else [v[1] for v in axes_main]

        if self.subset is False:
            cell_ann_axes = axes_main[-1][1:] if self.tree is None else axes_main[-1][2:]
            subset_ax = None
        else:
            cell_ann_axes = axes_main[-1][1:-1] if self.tree is None else axes_main[-1][2:-1]
            subset_ax = axes_main[-1][-1] if self.tree is None else axes_main[-1][-1]

        return tree_ax, heatmap_ax, cell_ann_axes, bin_ann_axes, subset_ax

    def setup_legend_axes(self, fig, disable_axis=False):
        num_cell_annotations = self.num_cell_annotations + 1 if self.subset else self.num_cell_annotations

        axes_legends = fig.subplots(
            nrows=1,
            ncols=1 + self.num_bin_annotations + num_cell_annotations,
            squeeze=False
        )[0]

        if disable_axis:
            for ax in axes_legends:
                ax.set_axis_off()
                ax.set_alpha(0.0)
                ax.patch.set_alpha(0.0)

        heatmap_ax_legend = axes_legends[0]
        cell_ann_axes_legends = None
        subset_ax_legend = None
        if self.cell_annotation_data is not None:
            if self.subset:
                cell_ann_axes_legends = axes_legends[1:num_cell_annotations]
                subset_ax_legend = axes_legends[num_cell_annotations]
            else:
                cell_ann_axes_legends = axes_legends[1:num_cell_annotations + 1]
                subset_ax_legend = None

        bin_ann_axes_legends = None
        if self.bin_annotation_data is not None:
            bin_ann_axes_legends = axes_legends[1 + num_cell_annotations:]

        return heatmap_ax_legend, cell_ann_axes_legends, bin_ann_axes_legends, subset_ax_legend

    def setup_primary_fig(self, fig):

        # split into main figure and legend
        fig_main, fig_legends = fig.subfigures(nrows=2, ncols=1, height_ratios=[5, 1], squeeze=True)
        fig_legends.patch.set_alpha(0.0)

        return fig_main, fig_legends

    def plot_annotation(self, data, title, bar_ax, legend_ax, horizontal=False, color_levels=None):

        if color_levels is None:
            im = bar_ax.imshow(data, aspect='auto', interpolation='none', cmap='Reds')
            legend_ax.grid(False)
            legend_ax.set_xticks([])
            legend_ax.set_yticks([])

            axins = legend_ax.inset_axes([0.5, 0.1, 0.05, 0.8])

            cbar = plt.colorbar(im, cax=axins)
            axins.set_title(title, fontsize='6')
            cbar.ax.tick_params(labelsize='4')
        else:
            bar_ax.imshow(data, aspect='auto', interpolation='none')
            levels = []
            patches = []
            for s, h in color_levels.items():
                levels.append(s)
                patches.append(Patch(facecolor=h, edgecolor=h))
            ncol = min(3, int(len(levels) ** (1 / 2)))
            legend_ax.legend(patches, levels, ncol=ncol,
                             frameon=True, loc=2, bbox_to_anchor=(0., 1.),
                             facecolor='white', edgecolor='white', fontsize='4',
                             title=title, title_fontsize='6')

        bar_ax.grid(False)

        if horizontal:
            bar_ax.set_yticks([0.], [title], rotation=0, fontsize='6')
            bar_ax.tick_params(axis='x', left=False, right=False)
        else:
            bar_ax.set_xticks([0.], [title], rotation=90, fontsize='6')
            bar_ax.tick_params(axis='y', left=False, right=False)

    def plot_cell_annotations(self):

        if self.cell_annotation_data is None:
            return

        for ax, ax_legend, annotation_field in zip(self.cell_ann_axes, self.cell_ann_axes_legends,
                                                   self.cell_annotation_data.columns):

            if self.cell_annotation_data[annotation_field].dtype.name in ('category', 'object', 'bool'):
                annotation_map = self.cell_annotation_data[annotation_field].to_dict()
                color_levels = colors.map_colormap_to_levels(set(annotation_map.values()))
                values = [[color_levels[annotation_map[cell]]] for cell in self.cell_order]
            else:
                values = np.array([self.cell_annotation_data[annotation_field].values])
                color_levels = None

            self.plot_annotation(
                values, annotation_field, ax, ax_legend, horizontal=False, color_levels=color_levels
            )

    def plot_bin_annotations(self):

        if self.bin_annotation_data is None:
            return

        for ax, ax_legend, annotation_field in zip(self.bin_ann_axes, self.bin_ann_axes_legends,
                                                   self.bin_annotation_data.columns):

            if self.bin_annotation_data[annotation_field].dtype.name in ('category', 'object', 'bool'):
                annotation_map = self.bin_annotation_data[annotation_field].to_dict()
                color_levels = colors.map_colormap_to_levels(
                    set(annotation_map.values()), colors=annotation_field
                )
                values = [color_levels[annotation_map[bin]] for bin in self.bins]
                values = np.array([values])
            else:
                values = np.array([self.bin_annotation_data[annotation_field].values])
                color_levels = None

            self.plot_annotation(
                values, annotation_field, ax, ax_legend, horizontal=True, color_levels=color_levels
            )

    def plot_phylogenetic_tree(self):
        if self.tree is None:
            return

        self.tree_ax.spines['top'].set_visible(False)
        self.tree_ax.spines['right'].set_visible(False)
        self.tree_ax.spines['bottom'].set_visible(True)
        self.tree_ax.spines['left'].set_visible(False)
        with plt.rc_context({'lines.linewidth': 0.5}):
            Bio.Phylo.draw(self.tree, label_func=lambda a: '', axes=self.tree_ax, do_show=False)
        self.tree_ax.tick_params(axis='x', labelsize=6)
        self.tree_ax.set_xlabel('branch length', fontsize=8)
        self.tree_ax.set_ylabel('')
        self.tree_ax.set_yticks([])

    def view_by_cell_annotation(self, annotation_label, annotation_value):

        cells_of_interest = self.cell_annotation_data[self.cell_annotation_data[annotation_label] == annotation_value]
        cells_of_interest = cells_of_interest.index
        indices = [self.cell_order.index(cell) for cell in cells_of_interest]

        sorted_indices = sorted(indices)
        self.add_subset_annotation(start=sorted_indices[0], end=sorted_indices[-1])

        current_lims = self.full_ylims
        inverted = False if current_lims[0] < current_lims[-1] else True
        current_lims = sorted(current_lims)
        current_lims = np.arange(current_lims[0], current_lims[-1], 1)

        start = current_lims[indices[0]]
        end = current_lims[indices[-1]]

        if inverted:
            start += 1
        else:
            end += 1

        self.heatmap_ax.set_ylim(start, end)
        for ax in self.cell_ann_axes:
            ax.set_ylim(start, end)
        self.subset_ax.set_ylim(start, end)

        return self.fig

    def show_cell_labels(self):
        current_lims = self.heatmap_ax.get_ylim()
        self.heatmap.show_cell_labels()
        self.heatmap_ax.set_ylim(*current_lims)
        return self.fig

    def hide_cell_labels(self):
        current_lims = self.heatmap_ax.get_ylim()
        self.heatmap.hide_cell_labels()
        self.heatmap_ax.set_ylim(*current_lims)
        return self.fig

    def reset_cells(self):
        self.add_subset_annotation()

        start, end = self.full_ylims
        self.heatmap_ax.set_ylim(start, end)
        for ax in self.cell_ann_axes:
            ax.set_ylim(start, end)
        self.subset_ax.set_ylim(start, end)

        return self.fig

    def add_subset_annotation(self, start=None, end=None):

        start = 0 if start is None else start
        end = len(self.cell_order) - 1 if end is None else end

        num_per_subset = abs(start - end) // 10

        subsets = {}
        curr_label = [0, 0]
        for i, cell in enumerate(self.cell_order):
            if i < start or i > end:
                subsets[cell] = -1
            else:
                subsets[cell] = curr_label[1]
                if curr_label[0] >= num_per_subset:
                    curr_label[0] = 0
                    curr_label[1] += 1
                else:
                    curr_label[0] += 1

        subset_data = self.cell_annotation_data.index.to_series().map(subsets)
        subset_data = subset_data.astype('category')

        self.cell_annotation_data['subset'] = subset_data

        annotation_map = subset_data.to_dict()
        color_levels = colors.map_colormap_to_levels(set(annotation_map.values()))
        values = [[color_levels[annotation_map[cell]]] for cell in self.cell_order]

        self.plot_annotation(
            values, 'subset', self.subset_ax, self.subset_ax_legend, horizontal=False, color_levels=color_levels
        )


def plot_cell_cn_matrix(
        adata: AnnData,
        layer_name='state',
        cell_order_fields=(),
        ax=None,
        max_cn=13,
        cmap=None,
        show_cell_ids=False):
    if len(cell_order_fields) > 0:
        cell_order_fields = reversed(list(cell_order_fields))
        cell_order_values = adata.obs[cell_order_fields].values.transpose()
        cell_ordering = adata.obs.index[np.lexsort(cell_order_values)]
    else:
        cell_ordering = None

    # could remove, current CN colormap can handle this
    X = adata.layers[layer_name]
    if max_cn is not None:
        X[X > max_cn] = max_cn

    df = pd.DataFrame(X, columns=adata.var.index, index=adata.obs.index)

    hmap = HeatMap(df, cell_order=cell_ordering, ax=ax, cmap=cmap)

    if show_cell_ids is True:
        hmap.show_cell_labels()
    else:
        hmap.hide_cell_labels()


def plot_cell_cn_matrix_fig(
        adata: AnnData,
        layer_name='state',
        tree=None,
        cell_order_fields=None,
        annotation_fields=None,
        var_annotation_fields=None,
        fig=None,
        cmap=None,
        max_cn=13,
        show_cell_ids=False,
        show_subsets=False
):
    if len(cell_order_fields) > 0:
        cell_order_fields = reversed(list(cell_order_fields))
        cell_order_values = adata.obs[cell_order_fields].values.transpose()
        cell_ordering = adata.obs.index[np.lexsort(cell_order_values)]
    else:
        cell_ordering = None

    # could remove, current CN colormap can handle this
    X = adata.layers[layer_name]
    if max_cn is not None:
        X[X > max_cn] = max_cn

    cell_annotation_data = adata.obs[annotation_fields]
    cell_annotation_data['cell_id'] = cell_annotation_data.index

    bin_annotation_data = adata.var[var_annotation_fields]
    bin_annotation_data['bin'] = bin_annotation_data.index

    df = pd.DataFrame(X, columns=adata.var.index, index=adata.obs.index)

    hmap = AnnotatedHeatMap(
        df,
        cell_order=cell_ordering,
        fig=fig,
        cmap=cmap,
        cell_annotation_data=cell_annotation_data,
        bin_annotation_data=bin_annotation_data,
        phylogenetic_tree=tree,
        subset=show_subsets
    )

    if show_cell_ids is True:
        hmap.show_cell_labels()
    else:
        hmap.hide_cell_labels()
