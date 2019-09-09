import logging
import yaml
import os
import pandas as pd
import numpy as np

import scgenome.utils
import scgenome.loaders.utils


def load_breakpoint_annotation_data(
        results_dir,
    ):
    """ Load breakpoint data from a pseudobulk run.

    Args:
        results_dir (str): results directory
    """
    bpcols = "prediction_id\tchromosome_1\tstrand_1\tposition_1\tchromosome_2\tstrand_2\t\
    position_2\thomology\tnum_split\tinserted\tmate_score\ttemplate_length_1\tlog_cdf\t\
    template_length_2\tlog_likelihood\ttemplate_length_min\tnum_reads\tnum_unique_reads\t\
    type\tnum_inserted\tsequence\tgene_id_1\tgene_name_1\tgene_location_1\tgene_id_2\t\
    gene_name_2\tgene_location_2\tdgv_ids\tis_germline\tis_dgv\tnum_patients\tis_filtered\t\
    dist_filtered\tbalanced\trearrangement_type".replace(' ', '').split('\t')

    suffix = 'destruct.csv.gz'
    if scgenome.loaders.utils.get_version(results_dir) == 'v0.2.11':
        suffix = 'destruct.tsv'
        bpcols = None

    breakpoint_data = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(results_dir, suffix):
        data = pd.read_csv(
            filepath, sep='\t', names=bpcols,
            dtype={
                'chromosome_1': 'str',
                'chromosome_2': 'str',
            })
        data['library_id'] = library_id
        data['sample_id'] = sample_id
        breakpoint_data.append(data)

    if len(breakpoint_data) == 0:
        return pd.DataFrame(), pd.DataFrame()

    breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

    return breakpoint_data


def load_breakpoint_count_data(
        results_dir
    ):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        results_dir (str): results directory
    """
    suffix = 'cell_counts_destruct.csv.gz'
    cols = ['cluster_id', 'cell_id', 'read_count']

    if scgenome.loaders.utils.get_version(results_dir) == 'v0.2.11':
        suffix = 'cell_counts_destruct.csv'
        cols = None

    breakpoint_count_data = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(results_dir, suffix):
        data = pd.read_csv(filepath, names=cols)
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
        results_dir (str): results directory
    """
    breakpoint_data = load_breakpoint_annotation_data(results_dir)
    breakpoint_count_data = load_breakpoint_count_data(results_dir)

    return {
        'breakpoint_data': breakpoint_data,
        'breakpoint_count_data': breakpoint_count_data,
    }


def load_cached_breakpoint_data(
        ticket_id,
        local_cache_directory,
    ):
    """ Load breakpoint tables from the cache
    
    Args:
        ticket_id (str): jira ticket for the analyis producing the results.
        local_cache_directory (str): local cache directory to search for results.
    
    KwArgs:
        additional_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    ticket_results_dirs = scgenome.loaders.utils.find_results_directories(
        ticket_id, local_cache_directory)

    if 'pseudobulk' not in ticket_results_dirs:
        raise ValueError(f'no pseudobulk found for ticket {ticket_id}')

    return load_breakpoint_data(
        ticket_results_dirs['pseudobulk'],
    )

