import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.sparse import issparse

from anndata import AnnData

import scgenome.cnplot


def plot_cn_profile(
        adata: AnnData,
        obs_id: str,
        value_layer_name=None,
        state_layer_name=None,
        ax=None,
        max_cn=13,
        chromosome=None,
        s=5,
        squashy=False,
        rawy=False,
    ):
    """Plot scatter points of copy number across the genome or a chromosome.

    Parameters
    ----------
    adata : AnnData
        copy number data
    obs_id : str
        observation to plot
    value_layer_name : str, optional
        layer with values for y axis, None for X, by default None
    state_layer_name : str, optional
        layer with states for colors, None for no color by state, by default None
    ax : [type], optional
        existing axis to plot into, by default None
    max_cn : int, optional
        max copy number for y axis, by default 13
    chromosome : [type], optional
        single chromosome plot, by default None
    s : int, optional
        size of scatter points, by default 5
    squashy : bool, optional
        compress y axis, by default False
    rawy : bool, optional
        raw data on y axis, by default False

    Examples
    -------

    .. plot::
        :context: close-figs

        import scgenome
        adata = scgenome.datasets.OV2295_HMMCopy_reduced()
        scgenome.pl.plot_cn_profile(adata, 'SA922-A90554B-R27-C43', value_layer_name='copy', state_layer_name='state')

    TODO: missing return
    """

    cn_data = adata.var.copy()

    if value_layer_name is not None:
        cn_value = adata[[obs_id], :].layers[value_layer_name]
    else:
        cn_value = adata[[obs_id], :].X

    if issparse(cn_value):
        cn_data['value'] = cn_value.toarray()[0]

    else:
        cn_data['value'] = np.array(cn_value)[0]

    # TODO: what if state is sparse
    cn_field_name = None
    if state_layer_name is not None:
        cn_data['state'] = np.array(adata[[obs_id], :].layers[state_layer_name][0])
        cn_field_name = 'state'

    if ax is None:
        ax = plt.gca()

    cn_data = cn_data.dropna(subset=['value'])

    scgenome.cnplot.plot_cell_cn_profile(
        ax, cn_data, 'value', cn_field_name=cn_field_name, max_cn=max_cn,
        chromosome=chromosome, s=s, squashy=squashy, rawy=rawy)

    return ax


def plot_var_profile(adata, value_field_name, cn_field_name=None, ax=None, max_cn=12, chromosome=None, s=5, squashy=False, rawy=False):
    """Plot scatter points of copy number across the genome or a chromosome.

    Parameters
    ----------
    adata : AnnData
        copy number data
    value_field_name : str, optional
        var field with values for y axis, None for X, by default None
    cn_field_name : str, optional
        var field with states for colors, None for no color by state, by default None
    ax : [type], optional
        existing axis to plot into, by default None
    max_cn : int, optional
        max copy number for y axis, by default 13
    chromosome : [type], optional
        single chromosome plot, by default None
    s : int, optional
        size of scatter points, by default 5
    squashy : bool, optional
        compress y axis, by default False
    rawy : bool, optional
        raw data on y axis, by default False

    Examples
    -------

    .. plot::
        :context: close-figs

        import scgenome
        adata = scgenome.datasets.OV2295_HMMCopy_reduced()
        scgenome.pl.plot_var_profile(adata[:, adata.var['gc'] > 0], 'gc', rawy=True)

    TODO: missing return
    """

    data = adata.var.copy()

    if ax is None:
        ax = plt.gca()

    scgenome.cnplot.plot_cell_cn_profile(
        ax, data, value_field_name, cn_field_name=cn_field_name, max_cn=max_cn,
        chromosome=chromosome, s=s, squashy=squashy, rawy=rawy)

    return ax

