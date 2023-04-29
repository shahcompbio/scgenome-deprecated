import seaborn as sns

import scgenome.tools.getters


def plot_gc_reads(adata, obs_id):
    """ Plot scatter points of gc by read count.

    Args:
        adata (anndata.AnnData): copy number
        obs_id (str): cell or clone to plot
    """

    data = scgenome.tools.getters.get_obs_data(adata, cell_id, var_columns=['gc'], layer_names=[None])
    data = data.rename(columns={'_X': 'reads'})
    data = data[data['gc'] > 0]

    sns.scatterplot(x='gc', y='reads', data=data, s=10)
