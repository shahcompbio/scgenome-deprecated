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

_table_suffixes_v0_0_0 = (
    ('align_metrics', '_alignment_metrics.csv.gz'),
)

_table_suffixes_v0_2_25 = (
    ('align_metrics', '_alignment_metrics.csv.gz'),
    ('gc_metrics', '_gc_metrics.csv.gz'),
)

table_suffixes = defaultdict(lambda: _table_suffixes_v0_2_25, {
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
})


def load_align_data_from_files(align_metrics, gc_metrics=None):
    results_tables = dict()

    results_tables["align_metrics"] = process_alignment_data(align_metrics, "align_metrics")

    if gc_metrics:
        results_tables["gc_metrics"] = process_alignment_data(gc_metrics, "gc_metrics")

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_align_data(
        results_dir,
):
    """ Load copy number tables
    
    Args:
        results_dir (str): results directory to load from.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'align' not in analysis_dirs:
        raise ValueError(f'no align found for directory {results_dir}')

    align_results_dir = analysis_dirs['align']
    assert len(align_results_dir) == 1
    align_results_dir = align_results_dir[0]

    manifest_filename = os.path.join(align_results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))

    # KLUDGE: 0.3.1 -> v0.3.1
    if not manifest['meta']['version'].startswith('v'):
        manifest['meta']['version'] = 'v' + manifest['meta']['version']

    version = manifest['meta']['version']

    results_tables = {}

    for table_name, suffix in table_suffixes[version]:
        filenames = scgenome.loaders.utils.find_filenames(manifest['filenames'], suffix)

        if len(filenames) != 1:
            raise ValueError(f'found filenames {filenames} for suffix {suffix}')

        filename = filenames[0]

        filepath = os.path.join(align_results_dir, filename)

        results_tables[table_name] = process_alignment_data(filepath, table_name)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


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
