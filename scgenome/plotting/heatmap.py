import Bio.Phylo
import colors
import numpy as np
import pandas as pd
import refgenome
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

cellids = list(df['cell_id'].unique())
annotations = pd.DataFrame(
    [(v, True if i < 3000 else False, False if i < 3000 else True, i) for i, v in enumerate(cellids)],
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
hmap = HeatMap(df, cell_annotation_data=annotations, bin_annotation_data=bin_annotations, cell_order=cellids[::-1])

plt.savefig('out.png')

"""


class HeatMap(object):
    def __init__(
            self,
            data,
            cell_annotation_data=None,
            bin_annotation_data=None,
            ax=None,
            cmap=None,
            cell_order=None,
            phylogenetic_tree=None
    ):

        self.sanity_check_inputs(
            cell_annotation_data, bin_annotation_data, ax, cell_order, phylogenetic_tree
        )

        self.data, self.cell_order, self.bins = self._reformat_data(data, phylogenetic_tree, cell_order)

        self.cmap = cmap
        self.tree = phylogenetic_tree

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

        if ax is not None:
            self.heatmap_ax = ax
            self.heatmap_ax.set_axis_on()
        else:
            main_fig, legend_fig = self.setup_primary_fig()
            self.tree_ax, self.heatmap_ax, self.cell_ann_axes, self.bin_ann_axes = self.setup_heatmap_axes(
                main_fig, disable_axis=True)
            self.heatmap_ax_legend, self.cell_ann_axes_legends, self.bin_ann_axes_legends = self.setup_legend_axes(
                legend_fig, disable_axis=True)

            self.heatmap_ax.set_axis_on()

            for cell_annotation in self.cell_ann_axes:
                cell_annotation.set_axis_on()

            for bin_annotation in self.bin_ann_axes:
                bin_annotation.set_axis_on()

        self.plot_heatmap()

        self.plot_cell_annotations()

        self.plot_bin_annotations()

        self.plot_phylogenetic_tree()

    def sanity_check_inputs(
            self, cell_annotation_data, bin_annotation_data,
            axes, cell_order, tree
    ):

        if axes is not None:
            assert cell_annotation_data is None, 'annotation bars not allowed when an axis object is supplied'
            assert bin_annotation_data is None, 'annotation bars not allowed when an axis object is supplied'
            assert tree is None, 'annotation bars not allowed when an axis object is supplied'

        if tree is not None and cell_order is not None:
            raise Exception('both tree and cell_order provided. please provide only one')

    def _reformat_data(self, data, tree, cell_order):
        data = data.fillna(0)

        if tree is not None and cell_order is not None:
            raise Exception('both tree and cell_order provided. please provide only one')

        if cell_order is None and tree is not None:
            cell_order = [a.name for a in self.tree.get_terminals()]
        elif cell_order is None:
            cell_order = list(self.data.index)

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

    def setup_heatmap_axes(self, fig, disable_axis=False):
        width_ratios = [1.0] + [0.02] * self.num_cell_annotations
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
        cell_ann_axes = axes_main[-1][1:] if self.tree is None else axes_main[-1][2:]
        bin_ann_axes = [v[0] for v in axes_main] if self.tree is None else [v[1] for v in axes_main]

        return tree_ax, heatmap_ax, cell_ann_axes, bin_ann_axes

    def setup_legend_axes(self, fig, disable_axis=False):

        axes_legends = fig.subplots(
            nrows=1,
            ncols=1 + self.num_bin_annotations + self.num_cell_annotations,
            squeeze=False
        )[0]

        if disable_axis:
            for ax in axes_legends:
                ax.set_axis_off()
                ax.set_alpha(0.0)
                ax.patch.set_alpha(0.0)

        heatmap_ax_legend = axes_legends[0]
        cell_ann_axes_legends = None
        if self.cell_annotation_data is not None:
            cell_ann_axes_legends = axes_legends[1:self.num_cell_annotations + 1]

        bin_ann_axes_legends = None
        if self.bin_annotation_data is not None:
            bin_ann_axes_legends = axes_legends[1 + self.num_cell_annotations:]

        return heatmap_ax_legend, cell_ann_axes_legends, bin_ann_axes_legends

    def setup_primary_fig(self):

        # split into main figure and legend
        fig = plt.figure(figsize=(10, 10))
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

                color_levels = colors.map_colormap_to_levels(self.cell_annotation_data[annotation_field].unique())

                annotation_map = self.cell_annotation_data[annotation_field].to_dict()

                color_values = [color_levels[annotation_map[cell]] for cell in self.cell_order]
                color_values = np.array([color_values])

                self.plot_annotation(
                    color_values, annotation_field, ax, ax_legend, horizontal=False, color_levels=color_levels
                )
            else:
                values = np.array([self.cell_annotation_data[annotation_field].values]).T
                self.plot_annotation(values, annotation_field, ax, ax_legend, horizontal=False)

    def plot_bin_annotations(self):

        if self.bin_annotation_data is None:
            return

        for ax, ax_legend, annotation_field in zip(self.bin_ann_axes, self.bin_ann_axes_legends,
                                                   self.bin_annotation_data.columns):

            if self.bin_annotation_data[annotation_field].dtype.name in ('category', 'object', 'bool'):

                color_category = annotation_field if annotation_field == 'cyto_band_giemsa_stain' else None

                color_levels = colors.map_colormap_to_levels(
                    self.bin_annotation_data[annotation_field].unique(), colors=color_category
                )

                annotation_map = self.bin_annotation_data[annotation_field].to_dict()

                color_values = [color_levels[annotation_map[bin]] for bin in self.bins]
                color_values = np.array([color_values])

                self.plot_annotation(
                    color_values, annotation_field, ax, ax_legend, horizontal=True, color_levels=color_levels
                )
            else:
                values = np.array([self.bin_annotation_data[annotation_field].values])
                self.plot_annotation(values, annotation_field, ax, ax_legend, horizontal=True)

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

    def plot_heatmap(self):

        mat = np.array(self.data)
        mat = np.nan_to_num(mat, nan=0)

        X_colors = colors.map_cn_colors_to_matrix(mat)

        self.heatmap_ax.imshow(X_colors, aspect='auto', interpolation='none')

        chridxs = np.array(self.bins.str.split(':').str[0])
        mat_chrom_idxs = np.array(sorted(np.unique(chridxs, return_inverse=True)[1]))

        chrom_boundaries = np.array(
            [0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [mat_chrom_idxs.shape[0] - 1]
        )

        chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
        chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2
        ordered_mat_chrom_idxs = mat_chrom_idxs[np.where(np.array([1] + list(np.diff(mat_chrom_idxs))) != 0)]
        chrom_names = np.array(refgenome.plot_chromosomes())[ordered_mat_chrom_idxs]

        self.heatmap_ax.set_xticks(chrom_mids)
        self.heatmap_ax.set_xticklabels(chrom_names, fontsize='6')
        self.heatmap_ax.set_xlabel('chromosome', fontsize=8)

        # cell labels
        self.heatmap_ax.set(yticks=range(len(self.data.index)))
        self.heatmap_ax.set(yticklabels=self.data.index.values)
        self.heatmap_ax.tick_params(axis='y', labelrotation=0)

        for val in chrom_boundaries[:-1]:
            self.heatmap_ax.axvline(x=val, linewidth=0.5, color='black', zorder=100)
