import seaborn as sns

import scgenome.tools.getters


def plot_gc_reads(adata, obs_id, **kwargs):
    """ Plot scatter points of gc by read count.

    Args:
        adata (anndata.AnnData): copy number
        obs_id (str): cell or clone to plot
    
    KwArgs:
        **kwargs: passed to sns.scatterplot
    """

    data = scgenome.tools.getters.get_obs_data(adata, obs_id, var_columns=['gc'], layer_names=[None])
    data = data.rename(columns={'_X': 'reads'})
    data = data[data['gc'] > 0]

    sns.scatterplot(x='gc', y='reads', data=data, **kwargs)
