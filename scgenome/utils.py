import pandas as pd
import numpy as np
import collections

from scgenome.constants import CN_DATA_ID, CELL_ID
from . import refgenome


chrom_names = refgenome.info.chromosomes

chrom_idxs = pd.Series(chrom_names)
chrom_idxs.name = 'chr'
chrom_idxs.index.name = 'chr_index'
chrom_idxs = chrom_idxs.reset_index()


def union_categories(dfs, cat_cols=None):
    """ Recreate specified categoricals on the union of categories inplace.

    Args:
        dfs (list of pandas.DataFrame): pandas dataframes to unify categoricals in-place.
    
    KwArgs:
        cat_cols (list of str): columns to unify categoricals, default None for any categorical.
    """

    # Infer all categorical columns if not given
    if cat_cols is None:
        cat_cols = set()
        for df in dfs:
            for col in df:
                if df[col].dtype.name == 'category':
                    cat_cols.add(col)

    # Get a list of categories for each column
    col_categories = collections.defaultdict(set)
    for col in cat_cols:
        for df in dfs:
            if col in df:
                col_categories[col].update(df[col].values)

    # Remove None and nan as they cannot be in the list of categories
    for col in col_categories:
        col_categories[col] = col_categories[col] - set([None, np.nan])

    # Create a pandas index for each set of categories
    for col, categories in col_categories.items():
        col_categories[col] = pd.Index(categories)

    # Set all categorical columns as having teh same set of categories
    for col in cat_cols:
        for df in dfs:
            if col in df:
                df[col] = df[col].astype('category')
                df[col] = df[col].cat.set_categories(col_categories[col])


def concat_with_categories(dfs, **kwargs):
    """ Concatenate dataframes retaining categorical columns
    """
    # Infer all categorical columns
    cat_cols = set()
    for df in dfs:
        for col in df:
            if df[col].dtype.name == 'category':
                cat_cols.add(col)

    # Get a list of categories for each column
    col_categories = collections.defaultdict(set)
    for df in dfs:
        for col in cat_cols:
            col_categories[col].update(df[col].values)

    # Remove None and nan as they cannot be in the list of categories
    for col in col_categories:
        col_categories[col] = col_categories[col] - set([None, np.nan])

    # Create a pandas index for each set of categories
    for col, categories in col_categories.items():
        col_categories[col] = pd.Index(categories)

    # Set all categorical columns as having teh same set of categories
    for df in dfs:
        for col in cat_cols:
            df[col] = df[col].astype('category')
            df[col] = df[col].cat.set_categories(col_categories[col])

    return pd.concat(dfs, **kwargs)


def cn_data_to_mat_data_ids(cn_data, data_id=CN_DATA_ID, cell_id=CELL_ID):
    matrix_data = (
        cn_data.set_index(['chr', 'start', 'cell_id'])
        [['reads', 'state', 'copy']]
            .unstack(level=2, fill_value=0.))
    copy = matrix_data[data_id].values
    measurement = copy.T

    cell_ids = matrix_data.columns.to_frame().loc[data_id][CELL_ID]
    return matrix_data, measurement, cell_ids


