import os
from collections import defaultdict

import pandas as pd
import scgenome.loaders.utils
import scgenome.utils
import yaml
from csverve.core import CsverveInput


_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]


def load_annotation_files(annotation_metrics):
    results_tables = {}

    results_tables['annotation_metrics'] = process_annotation_file(annotation_metrics)

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


def process_annotation_file(filepath):
    data = CsverveInput(filepath).read_csv()

    data['sample_id'] = [a.split('-')[-4] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[-3] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
