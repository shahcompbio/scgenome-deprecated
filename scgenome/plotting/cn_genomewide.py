"""

usage:
gwplot = GenomeWidePlot(
    data, col_to_plot, plot_func
)
plt.savefig(out_png)

data: dataframe with required cols: chr, start, end
col_to_plot: col in data frame with the data we want to plot
plot_func: function to be plotted


to access the axes object:
> gwplot.ax

(can also be provided when instantiating GenomeWidePlot)


To set x axis to one chromosome only:
> gwplot.set_axis_to_chromosome(chromosome)


Example usage:
data = csverve.read_csv("results/SHAH_H003347_T01_01_DLP01_hmmcopy_reads_small.csv.gz")

gwplot = GenomeWidePlot(
    data, 'copy', plt.scatter, color_col='state'
)

gwplot.set_axis_to_chromosome('1')
gwplot.squash_y_axis()
plt.savefig('out.png')



TODO : 'least surprise': would be better to have a supported list of plotting functions input as string: scatter, bar, line etc
and handle all of their peculiarities in this class instead of taking a callable func as input.
"""

import colors
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
from refgenome import chromosome_info
from refgenome import chromosomes
from refgenome import plot_chromosomes


class GenomeWidePlot(object):
    def __init__(
            self,
            data,
            y,
            ax=None,
            hue=None,
            kind='scatter',
            palette='cn'
    ):

        self.data = data
        self.y = y

        self.ax = plt.gca() if ax is None else ax

        assert kind in ('scatter', 'line', 'bar')
        self.kind = kind

        self.hue = hue

        self.cmap, self.color_reference = self.get_cmap(palette, hue)

        if isinstance(hue, str):
            self.c = self.data[hue]
        else:
            self.c = hue

        self.refgenome_chromosomes = chromosomes()
        self.refgenome_plot_chromosomes = plot_chromosomes()
        self.refgenome_chromosome_info = chromosome_info()

        self.add_chromosome_info()
        self.plot()
        self.setup_genome_axis()

    def get_cmap(self, palette, hue):
        matplotlib_cmaps = plt.colormaps()

        color_ref = None

        if palette == 'cn':
            palette = colors.Colors('cn').hex_color_reference

        if isinstance(palette, str):
            if palette in matplotlib_cmaps:
                palette = matplotlib.colors.Colormap(palette)
            else:
                palette = sns.color_palette(palette)
        elif isinstance(palette, dict):
            color_ref = palette

            if isinstance(hue, str):
                hue = self.data[hue].dropna().astype(int).values

            min_cn = int(hue.min())
            max_cn = int(hue.max())

            color_list = [palette[cn] for cn in range(min_cn, max_cn + 1)]
            palette = ListedColormap(color_list)
        else:
            raise NotImplementedError()

        return palette, color_ref

    def setup_genome_axis(self):
        self.ax.set_xlim((-0.5, self.refgenome_chromosome_info['chromosome_end'].max()))
        self.ax.set_xlabel('chromosome')
        self.ax.set_xticks([0] + list(self.refgenome_chromosome_info['chromosome_end'].values))
        self.ax.set_xticklabels([])
        self.ax.xaxis.tick_bottom()
        self.ax.yaxis.tick_left()
        self.ax.xaxis.set_minor_locator(
            matplotlib.ticker.FixedLocator(self.refgenome_chromosome_info['chromosome_mid'])
        )
        self.ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(self.refgenome_plot_chromosomes))
        sns.despine(ax=self.ax, offset=10, trim=True)

    def add_chromosome_info(self):
        self.data = self.data.merge(self.refgenome_chromosome_info)
        self.data = self.data[self.data['chr'].isin(chromosomes())]
        self.data['start'] = self.data['start'] + self.data['chromosome_start']
        self.data['end'] = self.data['end'] + self.data['chromosome_start']

    def set_y_ticks(self, max_cn=13):
        self.ax.set_ylim((-0.05 * max_cn, max_cn))
        self.ax.set_yticks(range(0, int(max_cn) + 1))
        self.ax.spines['left'].set_bounds(0, max_cn)

    def scatter_plot(self):

        if self.c is not None:
            plt.scatter(
                x=self.data['start'], y=self.data[self.y],
                c=self.c, cmap=self.cmap
            )
        else:
            plt.scatter(
                x=self.data['start'], y=self.data[self.y]
            )

        self.set_y_ticks()

    def line_plot(self):

        sns.lineplot(data=self.data, x='start', y=self.y, hue=self.hue, palette=self.cmap)

    def plot(self):
        if self.kind == 'scatter':
            self.scatter_plot()
        elif self.kind == 'line':
            sns.lineplot(data=self.data, x='start', y=self.y, hue=self.hue, palette=self.cmap)
        elif self.kind == 'bar':
            if self.hue:
                colors = [self.color_reference[v] for v in self.data[self.hue]]
                plt.bar(self.data['start'], self.data[self.y], color=colors)
            else:
                plt.bar(self.data['start'], self.data[self.y])
        else:
            raise NotImplementedError()

    def set_axis_to_chromosome(self, chromosome):
        chromosome_length = self.refgenome_chromosome_info.set_index('chr').loc[
            chromosome, 'chromosome_length']
        chromosome_start = self.refgenome_chromosome_info.set_index('chr').loc[chromosome, 'chromosome_start']
        chromosome_end = self.refgenome_chromosome_info.set_index('chr').loc[chromosome, 'chromosome_end']
        xticks = np.arange(0, chromosome_length, 2e7)
        xticklabels = ['{0:d}M'.format(int(x / 1e6)) for x in xticks]
        xminorticks = np.arange(0, chromosome_length, 1e6)
        self.ax.set_xlabel(f'chromosome {chromosome}')
        self.ax.set_xticks(xticks + chromosome_start)
        self.ax.set_xticklabels(xticklabels)
        self.ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xminorticks + chromosome_start))
        self.ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        self.ax.set_xlim((chromosome_start, chromosome_end))
        sns.despine(ax=self.ax, offset=10, trim=False)

    def squash_y_axis(self):

        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)

        self.ax.set_yscale('function', functions=(squash_f, squash_f))
