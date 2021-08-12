import os
from collections import defaultdict

import pandas as pd
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils
import yaml

_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]


def load_alignment_files(align_metrics, gc_metrics=None):
    results_tables = dict()

    results_tables["align_metrics"] = process_alignment_data(align_metrics, "align_metrics")

    if gc_metrics:
        results_tables["gc_metrics"] = process_alignment_data(gc_metrics, "gc_metrics")

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_alignment_results(results_dir):
    """ Load alignment metrics tables
    
    Args:
        results_dir (str): results directory to load from.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    alignment_metrics_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_alignment_metrics.csv.gz', analysis_type='alignment')

    gc_metrics_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_gc_metrics.csv.gz', analysis_type='alignment')

    return load_alignment_files(alignment_metrics_filepath, gc_metrics=gc_metrics_filepath)


def process_alignment_data(filepath, table_name):
    csv_input = scgenome.csvutils.CsvInput(filepath)

    dtypes_override = None
    if table_name == 'align_metrics':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'metrics_column_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}
    elif table_name == 'gc_metrics':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'alignment_gc_metrics_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}

    data = csv_input.read_csv(dtypes_override=dtypes_override)

    data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
