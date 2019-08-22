import pandas as pd
import numpy as np
import collections

from scgenome.constants import CN_DATA_ID, CELL_ID, VALUE_IDS, INDEX_IDS, \
    COPY_ID, NBIN_NCHR
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


def cn_data_to_mat_data_ids(cn_data, data_id=CN_DATA_ID, cell_id=CELL_ID,
                            index_ids=INDEX_IDS, value_ids=VALUE_IDS):
    matrix_data = (
        cn_data.set_index(index_ids)[value_ids]
            .unstack(level=2, fill_value=0.))
    copy = matrix_data[data_id].values
    measurement = copy.T

    cell_ids = matrix_data.columns.to_frame().loc[data_id][cell_id]
    return matrix_data, measurement, cell_ids


def cn_mat_to_cn_data(cn_mat, cell_id_vals=None, cell_id=CELL_ID,
                      value_id=COPY_ID):
    if cell_id_vals is None and cell_id not in cn_mat.columns:
        cell_id_vals = [f"cell{i}" for i in range(cn_mat.shape[0])]

    cn_mat[cell_id] = cell_id_vals
    hmm = pd.melt(cn_mat, id_vars=cell_id, var_name="chr_bin",
                  value_name=value_id)
    split = hmm["chr_bin"].str.split("_", expand=True)
    split.columns = ["chr", "bin"]
    split["bin"] = split["bin"].astype("int")
    split["bin"] = split["bin"].astype("int")
    hmm = hmm.drop(columns="chr_bin")
    hmm = pd.concat([split, hmm], axis=1)
    hmm["start"] = hmm["bin"] * 10
    hmm["end"] = hmm["start"] + 9
    hmm = hmm.sort_values([cell_id, "chr", "start"]).reset_index(drop=True)
    return hmm

# TODO this is used alot, implemented in different ways, might want combine
def cn_mat_as_df(cn_mat, chr_names):
    # TODO chr_names must be list of strings or we get an error
    num_sample = cn_mat.shape[0]
    num_bin = cn_mat.shape[1]

    if num_bin % len(chr_names) != 0:
        raise ValueError(NBIN_NCHR)
    bins_per_chr = num_bin / len(chr_names)

    chr = np.repeat(chr_names, bins_per_chr)
    ind = [str(b) for b in range(num_bin)]
    colnames = np.char.array(chr) + "_" + np.char.array(ind)
    rownames = [str(s) for s in range(num_sample)]
    return pd.DataFrame(data=cn_mat, columns=colnames, index=rownames)

