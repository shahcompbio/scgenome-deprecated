import os
from collections import defaultdict

import pandas as pd
import scgenome.loaders.utils
import scgenome.loaders.utils
import scgenome.utils
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

_table_suffixes_v0_0_0 = (
    ('hmmcopy_reads', '_multiplier0_reads.csv.gz'),
    ('hmmcopy_segs', '_multiplier0_segments.csv.gz'),
    ('hmmcopy_metrics', '_multiplier0_metrics.csv.gz'),
)

_table_suffixes_v0_2_25 = (
    ('hmmcopy_reads', '_reads.csv.gz'),
    ('hmmcopy_segs', '_segments.csv.gz'),
    ('hmmcopy_metrics', '_metrics.csv.gz'),
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
})

_type_fixes = [
    ('hmmcopy_reads', 'state'),
    ('hmmcopy_metrics', 'total_mapped_reads_hmmcopy'),
    ('hmmcopy_metrics', 'state_mode'),
]

def load_hmmcopy_data(
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

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'hmmcopy' not in analysis_dirs:
        raise ValueError(f'no hmmcopy found for directory {results_dir}')

    hmmcopy_results_dir = analysis_dirs['hmmcopy']

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
        filenames = scgenome.loaders.utils.find_filenames(manifest['filenames'], suffix)

        if len(filenames) != 1:
            raise ValueError(f'found filenames {filenames} for suffix {suffix}')

        filename = filenames[0]

        filepath = os.path.join(hmmcopy_results_dir, filename)

        usecols = None
        if table_name == 'hmmcopy_reads':
            usecols = hmmcopy_reads_cols

        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv(usecols=usecols)

        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

        for col in _categorical_cols:
            if col in data:
                data[col] = pd.Categorical(data[col])

        results_tables[table_name] = data

    # FIXUP: older hmmcopy results have total_mapped_reads instead of total_mapped_reads_hmmcopy
    results_tables['hmmcopy_metrics'] = results_tables['hmmcopy_metrics'].rename(
        columns={'total_mapped_reads': 'total_mapped_reads_hmmcopy'})

    # Some Hmmcopy columns should be int, fill with 0 for cells with 0 reads
    for table_name, column in _type_fixes:
        if results_tables[table_name][column].dtype.name == 'float64':
            results_tables[table_name][column] = results_tables[table_name][column].fillna(0).astype(int)

    scgenome.utils.union_categories(results_tables.values())

    _table_fixes[version](results_tables)

    return results_tables
