import os
from collections import defaultdict

import pandas as pd
import scgenome.loaders.utils
import scgenome.loaders.utils
import scgenome.utils
import scgenome.utils
import yaml

_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]

_table_suffixes_v0_0_0 = [
    ('align_metrics', '_alignment_metrics.csv.gz'),
]

_table_suffixes_v0_2_25 = [
    ('align_metrics', '_alignment_metrics.csv.gz'),
    ('gc_metrics', '_gc_metrics.csv.gz'),
]

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
}
)


def _table_fixes_v0_0_0(results_tables):
    pass  # TODO


def _table_fixes_v0_2_25(results_tables):
    pass  # TODO


_table_fixes = defaultdict(lambda: _table_fixes_v0_2_25, {
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
)


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

    manifest_filename = os.path.join(align_results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

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

        data = pd.read_csv(filepath)

        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

        for col in _categorical_cols:
            if col in data:
                data[col] = pd.Categorical(data[col])

        results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    _table_fixes[version](results_tables)

    return results_tables
