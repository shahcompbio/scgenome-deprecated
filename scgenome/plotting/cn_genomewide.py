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

from collections import defaultdict

import csverve.api as csverve
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
            plot_function,
            ax=None,
            color_col=None
    ):
        self.data = data
        self.y = y
        self.plot_function = plot_function
        self.color_col = color_col

        self.ax = plt.gca() if ax is None else ax

        self.cmap = self.get_cmap()

        self.refgenome_chromosomes = chromosomes()
        self.refgenome_plot_chromosomes = plot_chromosomes()
        self.refgenome_chromosome_info = chromosome_info()

        self.add_chromosome_info()
        self.plot()
        self.setup_genome_axis()

    @property
    def color_reference(self):

        color_reference = {
            0: '#3182BD', 1: '#9ECAE1', 2: '#CCCCCC', 3: '#FDCC8A',
            4: '#FC8D59', 5: '#E34A33', 6: '#B30000', 7: '#980043',
            8: '#DD1C77', 9: '#DF65B0', 10: '#C994C7', 11: '#D4B9DA'
        }

        color_reference = defaultdict(lambda: '#D4B9DA', color_reference)

        return color_reference

    def get_cmap(self):
        data = self.data[self.color_col].dropna().astype(int).values
        min_cn = int(data.min())
        max_cn = int(data.max())
        color_list = [self.color_reference[cn] for cn in range(min_cn, max_cn + 1)]
        cmap = ListedColormap(color_list)
        return cmap

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

    def plot(self):

        c = self.data[self.color_col].apply(lambda x: self.color_reference[x])

        if self.color_col is not None:
            self.plot_function(
                x=self.data['start'], y=self.data[self.y],
                c=c, cmap=self.cmap
            )
        else:
            self.plot_function(
                x=self.data['start'], y=self.data[self.y]
            )

        self.set_y_ticks()

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

