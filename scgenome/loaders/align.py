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


def load_alignment_files(align_metrics, gc_metrics=None):
    results_tables = dict()

    results_tables["align_metrics"] = process_alignment_data(align_metrics)

    if gc_metrics:
        results_tables["gc_metrics"] = process_alignment_data(gc_metrics)

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
        results_dir, '_alignment_metrics.csv.gz', 'alignment_metrics', analysis_type='alignment')

    gc_metrics_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, '_gc_metrics.csv.gz', 'alignment_gc_metrics', analysis_type='alignment')

    return load_alignment_files(alignment_metrics_filepath, gc_metrics=gc_metrics_filepath)


def process_alignment_data(filepath):
    data = CsverveInput(filepath).read_csv()

    data.query(f"cell_id != 'reference'", inplace=True)

    data['sample_id'] = [a.split('-')[-4] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[-3] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
