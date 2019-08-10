import logging
import yaml
import os
import pandas as pd

import scgenome.utils


standard_hmmcopy_reads_cols = [
    'chr',
    'start',
    'end',
    'cell_id',
    'gc',
    'reads',
    'copy',
    'state',
]


_categorical_cols = [
    'cell_id',
    'chr',
    'sample_id',
    'library_id',
]


_table_suffixes_v0_0_0 = {
    ('hmmcopy_reads', '_multiplier0_reads.csv.gz'),
    ('hmmcopy_segs', '_multiplier0_segments.csv.gz'),
    ('hmmcopy_metrics', '_multiplier0_metrics.csv.gz'),
}


_table_suffixes_v0_2_25 = {
    ('hmmcopy_reads', '_reads.csv.gz'),
    ('hmmcopy_segs', '_segments.csv.gz'),
    ('hmmcopy_metrics', '_metrics.csv.gz'),
}


table_suffixes = {
    'v0.0.0': _table_suffixes_v0_0_0,
    'v0.1.5': _table_suffixes_v0_0_0,
    'v0.2.2': _table_suffixes_v0_0_0,
    'v0.2.3': _table_suffixes_v0_0_0,
    'v0.2.6': _table_suffixes_v0_0_0,
    'v0.2.7': _table_suffixes_v0_0_0,
    'v0.2.9': _table_suffixes_v0_0_0,
    'v0.2.10': _table_suffixes_v0_0_0,
    'v0.2.11': _table_suffixes_v0_0_0,
    'v0.2.15': _table_suffixes_v0_0_0,
    'v0.2.19': _table_suffixes_v0_0_0,
    'v0.2.20': _table_suffixes_v0_0_0,
    'v0.2.25': _table_suffixes_v0_2_25,
    'v0.3.0': _table_suffixes_v0_2_25,
    'v0.3.1': _table_suffixes_v0_2_25,
}


def _table_fixes_v0_0_0(results_tables):
    pass # TODO


def _table_fixes_v0_2_25(results_tables):
    pass # TODO


_table_fixes = {
    'v0.0.0': _table_fixes_v0_0_0,
    'v0.1.5': _table_fixes_v0_0_0,
    'v0.2.2': _table_fixes_v0_0_0,
    'v0.2.3': _table_fixes_v0_0_0,
    'v0.2.6': _table_fixes_v0_0_0,
    'v0.2.7': _table_fixes_v0_0_0,
    'v0.2.9': _table_fixes_v0_0_0,
    'v0.2.10': _table_fixes_v0_0_0,
    'v0.2.11': _table_fixes_v0_0_0,
    'v0.2.15': _table_fixes_v0_0_0,
    'v0.2.19': _table_fixes_v0_0_0,
    'v0.2.20': _table_fixes_v0_0_0,
    'v0.2.25': _table_fixes_v0_2_25,
    'v0.3.0': _table_fixes_v0_2_25,
    'v0.3.1': _table_fixes_v0_2_25,
}


def _find_filenames(filenames, suffix):
    return [f for f in filenames if f.endswith(suffix)]


def load_hmmcopy_data(
        hmmcopy_results_dir,
        additional_reads_cols=None,
    ):
    """ Load copy number tables
    
    Args:
        hmmcopy_results_dir (str): path to hmmcopy results.

    KwArgs:
        additional_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    hmmcopy_reads_cols = standard_hmmcopy_reads_cols.copy()
    if additional_reads_cols is not None:
        hmmcopy_reads_cols.extend(additional_reads_cols)

    manifest_filename = os.path.join(hmmcopy_results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

    # KLUDGE: 0.3.1 -> v0.3.1
    if not manifest['meta']['version'].startswith('v'):
        manifest['meta']['version'] = 'v' + manifest['meta']['version']

    version = manifest['meta']['version']

    results_tables = {}

    for table_name, suffix in table_suffixes[version]:
        filenames = _find_filenames(manifest['filenames'], suffix)

        if len(filenames) != 1:
            raise ValueError(f'found filenames {filenames} for suffix {suffix}')

        filename = filenames[0]

        filepath = os.path.join(hmmcopy_results_dir, filename)

        usecols = None
        if table_name == 'hmmcopy_reads':
            usecols = hmmcopy_reads_cols

        data = pd.read_csv(filepath, usecols=usecols)

        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

        for col in _categorical_cols:
            if col in data:
                data[col] = pd.Categorical(data[col])

        results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    _table_fixes[version](results_tables)

    return results_tables


def load_cached_hmmcopy_data(
        ticket_id,
        local_cache_directory,
        additional_reads_cols=None,
    ):
    """ Load copy number tables from the cache
    
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

    if 'hmmcopy' not in ticket_results_dirs:
        raise ValueError(f'no hmmcopy found for ticket {ticket_id}')

    return load_hmmcopy_data(
        ticket_results_dirs['hmmcopy'],
        additional_reads_cols=additional_reads_cols,
    )

