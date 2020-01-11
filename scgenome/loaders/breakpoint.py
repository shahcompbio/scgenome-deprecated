import logging
import yaml
import os
import pandas as pd
import numpy as np

import scgenome.utils
import scgenome.loaders.utils
import scgenome.csvutils


def load_breakpoint_annotation_data(
        pseudobulk_dir,
    ):
    """ Load breakpoint data from a pseudobulk run.

    Args:
        pseudobulk_dir (str): results directory
    """
    suffix = 'destruct_breakpoints.csv.gz'

    breakpoint_data = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, suffix):
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()
        data['library_id'] = library_id
        data['sample_id'] = sample_id
        breakpoint_data.append(data)

    if len(breakpoint_data) == 0:
        return pd.DataFrame(), pd.DataFrame()

    breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

    return breakpoint_data


def load_breakpoint_count_data(
        pseudobulk_dir
    ):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        pseudobulk_dir (str): results directory
    """
    suffix = 'destruct_cell_counts.csv.gz'

    breakpoint_count_data = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, suffix):
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()
        data['library_id'] = library_id
        data['sample_id'] = sample_id
        breakpoint_count_data.append(data)

    breakpoint_count_data = pd.concat(breakpoint_count_data, ignore_index=True)
    breakpoint_count_data = breakpoint_count_data.rename(columns={'cluster_id': 'prediction_id'})

    # KLUDGE: normal reads are not filtered properly, filter by their prefix
    breakpoint_count_data = breakpoint_count_data[~breakpoint_count_data['cell_id'].str.startswith('HS')]

    return breakpoint_count_data


def load_breakpoint_data(results_dir):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        results_dir (str): results directory to load from.
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['pseudobulk']

    elif 'breakpoint_calling' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['breakpoint_calling']

    else:
        raise ValueError(f'no breakpoints found for directory {results_dir}')

    breakpoint_data = load_breakpoint_annotation_data(pseudobulk_dir)
    breakpoint_count_data = load_breakpoint_count_data(pseudobulk_dir)

    # TODO: fix upstream
    for col in ('prediction_id', 'position_1', 'position_2', 'read_count'):
        for df in (breakpoint_data, breakpoint_count_data):
            if col in df:
                df[col] = df[col].astype(int)

    return {
        'breakpoint_data': breakpoint_data,
        'breakpoint_count_data': breakpoint_count_data,
    }


