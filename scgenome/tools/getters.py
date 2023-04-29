from anndata import AnnData
from pandas import DataFrame


def get_obs_data(
        adata: AnnData,
        obs_id: str,
        var_columns=None,
        layer_names=None) -> DataFrame:
    """ Get a dataframe for an observation including .var and specified layers

    Parameters
    ----------
    adata : AnnData
        anndata from which to retrieve observation data.
    obs_id : str
        Observation ID
    var_columns : list, optional
        List of var columns to add to returned dataframe, by default None, all columns
    layer_names : list, optional
        List of layer names to add to returned dataframe, by default None, all layers

    Returns
    -------
    DataFrame
        dataframe for an observation
    """
    data = adata.var

    if var_columns is not None:
        data = data[var_columns]

    if layer_names is None:
        layer_names = adata.layers.keys()

    for layer_name in layer_names:
        assert layer_name not in data, f'layer {layer_name} also in .var'
        layer_data = adata[obs_id].to_df(layer_name).T
        assert len(layer_data.columns) == 1
        if layer_name is not None:
            layer_data.columns = [layer_name]
        else:
            layer_data.columns = ['_X']
        data = data.merge(layer_data, left_index=True, right_index=True, how='left')
    
    return data

