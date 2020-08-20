import logging
import yaml
import os
import pandas as pd
import numpy as np

import scgenome.utils
import scgenome.loaders.utils
import scgenome.csvutils


def load_breakpoint_annotation_data(files, is_lumpy=False):
    """ Load breakpoint data from a pseudobulk run.

    Args:
        pseudobulk_dir (str): results directory
        suffix (str): suffix of breakpoint annotation tables
    """

    chrom_1_colname = "chromosome_1"
    chrom_2_colname = "chromosome_2"
    
    if is_lumpy:
        chrom_1_colname = "chrom1"
        chrom_2_colname = "chrom2"      

    breakpoint_data = []

    for sample_id, library_id, filepath in files:
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()

        data[chrom_1_colname] = data[chrom_1_colname].astype(str)
        data[chrom_2_colname] = data[chrom_2_colname].astype(str)

        if library_id is not None:
            data['library_id'] = library_id

        if sample_id is not None:
            data['sample_id'] = sample_id

        breakpoint_data.append(data)

    if len(breakpoint_data) == 0:
        return pd.DataFrame(), pd.DataFrame()

    breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

    return breakpoint_data


def load_breakpoint_count_data(files):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        pseudobulk_dir (str): results directory
        suffix (str): suffix of breakpoint count tables
    """
    breakpoint_count_data = []

    for sample_id, library_id, filepath in files:
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()

        if library_id is not None:
            data['library_id'] = pd.Series([library_id], dtype="category")

        if sample_id is not None:
            data['sample_id'] = pd.Series([sample_id], dtype="category")

        breakpoint_count_data.append(data)

    breakpoint_count_data = pd.concat(breakpoint_count_data, ignore_index=True)
    breakpoint_count_data = breakpoint_count_data.rename(columns={'cluster_id': 'prediction_id'})

    # KLUDGE: normal reads are not filtered properly, filter by their prefix, and having '-' in cell id
    breakpoint_count_data = breakpoint_count_data.loc[~breakpoint_count_data['cell_id'].str.startswith('HS'), :]
    breakpoint_count_data = breakpoint_count_data.loc[breakpoint_count_data['cell_id'].apply(lambda a: '-' in a), :]

    return breakpoint_count_data


def load_breakpoint_data_from_files(annotation_file, counts_file, lumpy=False):

    annotation_files  = scgenome.loaders.utils._prep_filenames_for_loading(annotation_file)

    breakpoint_data = load_breakpoint_annotation_data(annotation_files, is_lumpy=lumpy)

    count_files = scgenome.loaders.utils._prep_filenames_for_loading(counts_file)

    breakpoint_count_data = load_breakpoint_count_data(count_files)

    return _process_breakpoint_data(breakpoint_data, breakpoint_count_data)


def load_breakpoint_data(
        results_dir,
    ):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        results_dir (str): results directory to load from.
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        breakpoint_calling_dir = analysis_dirs['pseudobulk']
        annotation_suffix = 'destruct.csv.gz'
        count_suffix = 'cell_counts_destruct.csv.gz'

    elif 'breakpoint_calling' in analysis_dirs:
        breakpoint_calling_dir = analysis_dirs['breakpoint_calling']

        if len(breakpoint_calling_dir) == 0:
            raise ValueError(f'found {len(breakpoint_calling_dir)} dirs for breakpoint_calling')

        elif len(breakpoint_calling_dir) > 1:
            if filter_sample_id is None:
                raise ValueError(f'found {len(breakpoint_calling_dir)} without filter_sample_id')

            filtered_breakpoint_calling_dir = list(filter(lambda a: f'sample_{filter_sample_id}' in a, breakpoint_calling_dir))

            if len(filtered_breakpoint_calling_dir) != 1:
                raise ValueError(f'found {len(filtered_breakpoint_calling_dir)} in {breakpoint_calling_dir} matching filter_sample_id')

            breakpoint_calling_dir = filtered_breakpoint_calling_dir

        breakpoint_calling_dir = breakpoint_calling_dir[0]
        annotation_suffix = 'destruct_breakpoints.csv.gz'
        count_suffix = 'destruct_cell_counts.csv.gz'

    else:
        raise ValueError(f'no breakpoints found for directory {results_dir}')

        
    annotation_files  = scgenome.loaders.utils.get_pseudobulk_files(
        breakpoint_calling_dir, annotation_suffix)

    breakpoint_data = load_breakpoint_annotation_data(annotation_files)

    count_files = scgenome.loaders.utils.get_pseudobulk_files(
        breakpoint_calling_dir, count_suffix)

    breakpoint_count_data = load_breakpoint_count_data(count_files)

    return _process_breakpoint_data(breakpoint_data, breakpoint_count_data)


def _process_breakpoint_data(breakpoint_data, breakpoint_count_data):
    # TODO: fix upstream
    for col in ('prediction_id', 'position_1', 'position_2', 'read_count'):
        for df in (breakpoint_data, breakpoint_count_data):
            if col in df:
                df[col] = df[col].astype(int)

    return {
        'breakpoint_data': breakpoint_data,
        'breakpoint_count_data': breakpoint_count_data,
    }


