import click
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from scgenome import cncluster


@click.command()
@click.argument('cn_path', type=click.Path(exists=True))
@click.argument('cl_path', type=click.Path())
@click.option(
    '--plot', '-p', type=click.Path(), default=None,
    help='output plot file'
)
def main(cn_path, cl_path, plot):
    """
    Reduces dimensionality of cell copy number found in CN_PATH, cluster the
    cells using their lower dimension embedding, output table into CL_PATH.
    """

    cn = pd.read_table(
        cn_path, index_col=['chr', 'start', 'end', 'width'], dtype={'chr': str}
    )
    cl = cncluster.umap_hdbscan_cluster(cn)
    cl.to_csv(cl_path, sep='\t', index=False)

    if plot:
        fig, ax = plt.subplots(tight_layout=True)
        cncluster.plot_umap_clusters(ax, cl)
        plt.savefig(plot)


if __name__ == '__main__':
    main()
