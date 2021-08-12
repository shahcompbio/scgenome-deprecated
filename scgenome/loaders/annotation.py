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


def load_annotation_files(annotation_metrics):
    results_tables = {}

    results_tables['annotation_metrics'] = process_annotation_file(annotation_metrics, 'annotation_metrics')

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_annotation_results(
        results_dir,
):
    """ Load annotation tables
    
    Args:
        results_dir (str): results directory to load from.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    annotation_metrics_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_metrics.csv.gz', analysis_type='annotation')

    return load_annotation_files(annotation_metrics_filepath)


def process_annotation_file(filepath, table_name):
    csv_input = scgenome.csvutils.CsvInput(filepath)

    dtypes_override = None
    if table_name == 'annotation_metrics':
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'metrics_column_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}

    data = csv_input.read_csv(dtypes_override=dtypes_override)

    data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
