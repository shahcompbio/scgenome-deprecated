import pandas as pd
import numpy as np
import collections
from itertools import product

from scgenome.constants import CN_DATA_ID, CELL_ID, VALUE_IDS, INDEX_IDS, \
    COPY_ID, NBIN_NCHR, BHC_ID, ORIGIN_ID
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


def cn_data_to_mat_data_ids(cn_data, data_id=BHC_ID, cell_id=CELL_ID,
                            index_ids=INDEX_IDS, value_ids=VALUE_IDS):
    matrix_data = (
        cn_data.set_index(index_ids)[value_ids]
            .unstack(level=2, fill_value=0.))
    copy = matrix_data[data_id].values
    measurement = copy.T

    cell_ids = matrix_data.columns.to_frame().loc[data_id][cell_id]
    return matrix_data, measurement, cell_ids


def cn_data_to_mat1(cn_data, cluster_field_name="bhc_cluster_id",
                    origin_field_name=None):
    matrix_data = cn_data.merge(chrom_idxs)
    columns = ['chr_index', 'start', 'cell_id', cluster_field_name]
    levels = ['cell_id', cluster_field_name]
    if origin_field_name is not None:
        columns.append(origin_field_name)
        levels.append(origin_field_name)
    matrix_data = (matrix_data.set_index(columns)["state"]
                        .unstack(level=levels).fillna(0))

    copy = matrix_data.values
    measurement = copy.T

    cell_ids = matrix_data.columns.to_frame().loc["state"]["cell_id"]
    return matrix_data, measurement, cell_ids


def plot_get_mat(cn_data, cluster_field_name="bhc_cluster_id", origin_field_name=None,
                 cn_field_name="state"):
    plot_data = cn_data.merge(chrom_idxs)
    columns = ['chr_index', 'start', 'cell_id']
    levels = ['cell_id']
    if cluster_field_name is not None:
        columns.append(cluster_field_name)
        levels.append(cluster_field_name)
    if origin_field_name is not None:
        columns.append(origin_field_name)
        levels.append(origin_field_name)
    plot_data = (plot_data.set_index(columns)[cn_field_name]
                 .unstack(level=levels).fillna(0))
    return plot_data


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


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


def get_cn_data_submixture(cn_data, num_cells, hmmcopy_tickets,
                           proportions=None, origin_field="origin_id",
                           id_field_name="cell_id", seed=None):
    """

    :param hmmcopy_tickets: list of hmmcopy tickets
    :param sample_ids: list of lists of sample ids. first element corresponds
    to first element in `hmmcopy_tickets` etc.
    :return:
    """
    n_sample = len(hmmcopy_tickets)
    if proportions is None:
        proportions = [1/n_sample for i in range(n_sample)]

    cell_counts = (num_cells * np.array(proportions)).astype("int")

    sub_datasets = []
    for i in range(n_sample):
        jira_cn_data = cn_data[cn_data[origin_field] == hmmcopy_tickets[i]]
        sub_cn_data = subsample_cn_data(jira_cn_data, cell_counts[i],
                                        id_field_name, seed)
        sub_datasets.append(sub_cn_data)

    mixed = pd.concat(sub_datasets)
    return mixed


def subsample_cn_data(cn_data, num_cells, id_field_name="cell_id", seed=None):
    keep_ids = pd.Series(
        cn_data[id_field_name].unique()).sample(num_cells, random_state=seed)

    return cn_data[cn_data[id_field_name].isin(keep_ids)]


def get_mixture_labels(cn_data, obs_name="bhc_cluster_id",
                       exp_name="origin_id_int"):
    sub_cn_data = cn_data[["cell_id", obs_name, exp_name]].drop_duplicates()
    return sub_cn_data