import os
from collections import defaultdict

import pandas as pd
import scgenome.loaders.utils
import scgenome.utils
import scgenome.csvutils
import yaml


_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]


_table_suffixes_v0_2_25 = (
    ('annotation_metrics', '_metrics.csv.gz'),
)


table_suffixes = defaultdict(lambda: _table_suffixes_v0_2_25, {
    'v0.2.25': _table_suffixes_v0_2_25,
    'v0.3.0': _table_suffixes_v0_2_25,
    'v0.3.1': _table_suffixes_v0_2_25,
})


_dtype_fixes = [
#    ('annotation_metrics', 'total_mapped_reads_hmmcopy', 'Int64'),
#    ('annotation_metrics', 'is_s_phase_prob', 'float64'),
#    ('annotation_metrics', 'mad_neutral_state', 'float64'),
]


def load_annotation_data(
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

    if 'annotation' not in analysis_dirs:
        raise ValueError(f'no annotation found for directory {results_dir}')

    annotation_results_dir = analysis_dirs['annotation']

    manifest_filename = os.path.join(annotation_results_dir, 'metadata.yaml')
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

        filepath = os.path.join(annotation_results_dir, filename)

        csv_input = scgenome.csvutils.CsvInput(filepath)

        dtypes_override = None
        if table_name == 'annotation_metrics':
            dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
            dtypes_filename = os.path.join(dtypes_directory, 'metrics_column_defs.yaml')
            dtypes_override = yaml.load(open(dtypes_filename))
            dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}

        data = csv_input.read_csv(dtypes_override=dtypes_override)

        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

        for col in _categorical_cols:
            if col in data:
                data[col] = pd.Categorical(data[col])

        results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    for table_name, column_name, dtype in _dtype_fixes:
        if dtype == 'Int64':
            results_tables[table_name][column_name] = results_tables[table_name][column_name].astype('float').astype(dtype)
        else:
            results_tables[table_name][column_name] = results_tables[table_name][column_name].astype(dtype)

    if 'is_s_phase' in results_tables['annotation_metrics']:
        results_tables['annotation_metrics']['is_s_phase'] = results_tables['annotation_metrics']['is_s_phase'].fillna(False).astype(bool)

    return results_tables

