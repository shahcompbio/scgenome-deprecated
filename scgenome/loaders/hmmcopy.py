import os
from collections import defaultdict

import pandas as pd
import scgenome.loaders.utils
import scgenome.utils
import yaml
from csverve.core import CsverveInput


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


def load_hmmcopy_files(
        hmmcopy_reads, hmmcopy_segs, hmmcopy_metrics,
        additional_reads_cols=None
    ):
    results_tables = {}

    hmmcopy_reads_cols = standard_hmmcopy_reads_cols.copy()
    if additional_reads_cols is not None:
        hmmcopy_reads_cols.extend(additional_reads_cols)

    results_tables["hmmcopy_reads"] = process_hmmcopy_data(hmmcopy_reads, usecols=hmmcopy_reads_cols)
    results_tables["hmmcopy_segs"] = process_hmmcopy_data(hmmcopy_segs)
    results_tables["hmmcopy_metrics"] = process_hmmcopy_data(hmmcopy_metrics)

    # FIXUP: older hmmcopy results have total_mapped_reads instead of total_mapped_reads_hmmcopy
    results_tables['hmmcopy_metrics'] = results_tables['hmmcopy_metrics'].rename(
        columns={'total_mapped_reads': 'total_mapped_reads_hmmcopy'})

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_hmmcopy_results(
        results_dir,
        additional_reads_cols=None,
    ):
    """ Load copy number tables
    
    Args:
        results_dir (str): results directory to load from.

    KwArgs:
        additional_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    hmmcopy_reads_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_reads.csv.gz', analysis_type='hmmcopy')

    hmmcopy_segs_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_segments.csv.gz', analysis_type='hmmcopy')

    hmmcopy_metrics_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_metrics.csv.gz', analysis_type='hmmcopy')

    return load_hmmcopy_files(
        hmmcopy_reads_filepath, hmmcopy_segs_filepath, hmmcopy_metrics_filepath,
        additional_reads_cols=additional_reads_cols,
    )


def process_hmmcopy_data(filepath, usecols=None):
    data = CsverveInput(filepath).read_csv()

    data['sample_id'] = [a.split('-')[-4] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[-3] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
