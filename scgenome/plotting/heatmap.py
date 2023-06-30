import warnings

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

"""


class HeatMap(object):
    def __init__(
            self,
            data,
            cell_annotation_data=None,
            bin_annotation_data=None,
            ax=None,
            cmap=None,
            max_cn=None,
            cell_order=None,
            tree=None
    ):

        self.cmap = 'viridis' if cmap is None else cmap

        self.data = data
        self.data = self.data.fillna(0)
        if max_cn:
            self.data = self.data.clip(upper=max_cn)

        self.tree = tree

        if tree is not None:
            if cell_order is not None:
                warnings.warn('both tree and cell_order provided. tree takes precedence')

            self.cell_order = [a.name for a in tree.get_terminals()]

        elif cell_order is not None:
            self.cell_order = cell_order
        else:
            self.cell_order = list(self.data.index)

        self.data.index = pd.Categorical(self.data.index, self.cell_order)
        self.data = self.data.sort_index()

        self.bins = self.get_bins()
        self.sort_data()

        self.cell_annotation_data = cell_annotation_data
        if self.cell_annotation_data is not None:
            self.cell_annotation_data = self.cell_annotation_data.set_index('cell_id')

        self.bin_annotation_data = bin_annotation_data
        if self.bin_annotation_data is not None:
            self.bin_annotation_data = self.bin_annotation_data.set_index('bin')

        if ax is None:
            self.setup_axes()
        else:
            assert cell_annotation_data is None, 'annotation bars not allowed when an axis object is supplied'
            assert bin_annotation_data is None, 'annotation bars not allowed when an axis object is supplied'
            assert tree is None, 'annotation bars not allowed when an axis object is supplied'
            self.heatmap_ax = ax

    def setup_axes(self):

        num_phylo = 1 if self.tree is not None else 0
        heatmap_ax_col_idx = 1 if self.tree is not None else 0
        tree_ax_idx = 1 if self.tree is not None else 0
        num_annotations = len(self.cell_annotation_data.columns) if self.cell_annotation_data is not None else 0
        num_var_annotations = len(self.bin_annotation_data.columns) if self.bin_annotation_data is not None else 0
        heatmap_ax_row_idx = num_var_annotations + 1

        fig = plt.figure()
        fig_main, fig_legends = fig.subfigures(nrows=2, ncols=1, height_ratios=[5, 1], squeeze=True)
        fig_legends.patch.set_alpha(0.0)

        width_ratios = [0.5] * num_phylo + [1] + [0.005] + [0.02] * num_annotations
        height_ratios = [0.02] * num_var_annotations + [0.01] + [1]

        axes = fig_main.subplots(
            nrows=len(height_ratios), ncols=len(width_ratios),
            width_ratios=width_ratios, height_ratios=height_ratios,
            squeeze=False, gridspec_kw=dict(hspace=0.02, wspace=0.02))

        # Turn off axes for all annotation rows and columns
        for ax in axes[:1, :].flatten():
            ax.set_axis_off()
        for ax in axes[:, heatmap_ax_col_idx + 1:].flatten():
            ax.set_axis_off()

        # Re-enable axes and remove ticks for row annotations
        for ax in axes[heatmap_ax_row_idx, heatmap_ax_col_idx + 2:].flatten():
            ax.set_axis_on()
            ax.set_yticks([])

        # Re-enable axes and remove ticks for column annotations
        for ax in axes[:heatmap_ax_row_idx - 1, heatmap_ax_col_idx].flatten():
            ax.set_axis_on()
            ax.set_xticks([])

        axes_legends = fig_legends.subplots(
            nrows=1, ncols=1 + num_annotations + num_var_annotations, squeeze=False)[0]
        for ax in axes_legends:
            ax.set_axis_off()
            ax.set_alpha(0.0)
            ax.patch.set_alpha(0.0)

        self.tree_ax = axes[heatmap_ax_row_idx, tree_ax_idx]

        self.heatmap_ax = axes[heatmap_ax_row_idx, heatmap_ax_col_idx]
        self.heatmap_ax_legend = axes_legends[0]

        self.cell_annotation_axes = axes[heatmap_ax_row_idx, heatmap_ax_col_idx + 2:]
        self.cell_annotation_axes_legends = axes_legends[1:]

        self.bin_annotation_axes = axes[:, heatmap_ax_col_idx]
        self.bin_annotation_axes_legends = axes_legends[1 + num_annotations:]

    def get_bins(self):

        bins = self.data.columns

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

    def sort_data(self):
        self.data = self.data[self.bins]

    def plot_categorical_annotation(self, data, color_levels, title, bar_ax, legend_ax, horizontal=False):

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

    def plot_continuous_annotation(self, data, title, bar_ax, legend_ax, horizontal=False):
        im = bar_ax.imshow(data, aspect='auto', interpolation='none', cmap='Reds')

        bar_ax.grid(False)
        if horizontal:
            bar_ax.set_yticks([0.], [title], rotation=0, fontsize='6')
            bar_ax.tick_params(axis='x', left=False, right=False)
        else:
            bar_ax.set_xticks([0.], [title], rotation=90, fontsize='6')
            bar_ax.tick_params(axis='y', left=False, right=False)

        legend_ax.grid(False)
        legend_ax.set_xticks([])
        legend_ax.set_yticks([])

        axins = legend_ax.inset_axes([0.5, 0.1, 0.05, 0.8])

        cbar = plt.colorbar(im, cax=axins)
        axins.set_title(title, fontsize='6')
        cbar.ax.tick_params(labelsize='4')

    def plot_cell_annotations(self):

        if self.cell_annotation_data is None:
            return

        for ax, ax_legend, annotation_field in zip(self.cell_annotation_axes, self.cell_annotation_axes_legends,
                                                   self.cell_annotation_data.columns):

            if self.cell_annotation_data[annotation_field].dtype.name in ('category', 'object', 'bool'):

                color_levels = colors.map_colormap_to_levels(self.cell_annotation_data[annotation_field].unique())

                annotation_map = self.cell_annotation_data[annotation_field].to_dict()

                color_values = [color_levels[annotation_map[cell]] for cell in self.cell_order]
                color_values = np.array([color_values])

                self.plot_categorical_annotation(
                    color_values, color_levels, annotation_field, ax, ax_legend, horizontal=False
                )
            else:
                values = np.array([self.cell_annotation_data[annotation_field].values]).T
                self.plot_continuous_annotation(values, annotation_field, ax, ax_legend, horizontal=False)

    def plot_bin_annotations(self):

        if self.bin_annotation_data is None:
            return

        for ax, ax_legend, annotation_field in zip(self.bin_annotation_axes, self.bin_annotation_axes_legends,
                                                   self.bin_annotation_data.columns):

            if self.bin_annotation_data[annotation_field].dtype.name in ('category', 'object', 'bool'):

                color_category = annotation_field if annotation_field == 'cyto_band_giemsa_stain' else None

                color_levels = colors.map_colormap_to_levels(
                    self.bin_annotation_data[annotation_field].unique(), colors=color_category
                )

                annotation_map = self.bin_annotation_data[annotation_field].to_dict()

                color_values = [color_levels[annotation_map[bin]] for bin in self.bins]
                color_values = np.array([color_values])

                self.plot_categorical_annotation(
                    color_values, color_levels, annotation_field, ax, ax_legend, horizontal=True
                )
            else:
                values = np.array([self.bin_annotation_data[annotation_field].values])
                self.plot_continuous_annotation(values, annotation_field, ax, ax_legend, horizontal=True)

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

        # raise Exception(self.heatmap_ax)

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

    def main(self):

        self.plot_heatmap()

        self.plot_cell_annotations()

        self.plot_bin_annotations()

        self.plot_phylogenetic_tree()
