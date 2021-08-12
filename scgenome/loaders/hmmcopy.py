import os
from collections import defaultdict

import pandas as pd
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils
import yaml

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

    results_tables["hmmcopy_reads"] = process_hmmcopy_data(hmmcopy_reads, "hmmcopy_reads", usecols=hmmcopy_reads_cols)
    results_tables["hmmcopy_segs"] = process_hmmcopy_data(hmmcopy_segs, "hmmcopy_segs")
    results_tables["hmmcopy_metrics"] = process_hmmcopy_data(hmmcopy_metrics, "hmmcopy_metrics")

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


def process_hmmcopy_data(filepath, table_name, usecols=None):
    csv_input = scgenome.csvutils.CsvInput(filepath)

    dtypes_override = None
    if table_name == 'hmmcopy_metrics':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'metrics_column_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}
    elif table_name == 'hmmcopy_reads':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'hmmcopy_reads_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}
    elif table_name == 'hmmcopy_segs':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'hmmcopy_segments_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}

    data = csv_input.read_csv(usecols=usecols, dtypes_override=dtypes_override)

    data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
