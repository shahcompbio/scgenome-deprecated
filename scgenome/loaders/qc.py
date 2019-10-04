import os
import yaml
import pandas as pd

import scgenome.loaders.utils
from scgenome.loaders.align import load_align_data
from scgenome.loaders.hmmcopy import load_hmmcopy_data
from scgenome.loaders.annotation import load_annotation_data


def _calculate_annotation_metrics(results_tables):
    common_columns = list(set(results_tables['align_metrics'].columns).intersection(
        results_tables['hmmcopy_metrics'].columns))

    align_idx = results_tables['align_metrics'].set_index(common_columns).index
    hmmcopy_idx = results_tables['hmmcopy_metrics'].set_index(common_columns).index

    for idx in align_idx.difference(hmmcopy_idx):
        raise ValueError(f'found {idx} in align but not hmmcopy')

    for idx in hmmcopy_idx.difference(align_idx):
        raise ValueError(f'found {idx} in hmmcopy but not align')

    data = pd.merge(
        results_tables['align_metrics'],
        results_tables['hmmcopy_metrics'],
        on=common_columns,
        how='left',
    )

    return data


def load_cell_state_prediction(results_dir):
    """ Load cell state prediction table.

    Args:
        results_dir (str): results directory to load from.
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'cell_state_prediction' not in analysis_dirs:
        raise ValueError(f'no cell state prediction found for directory {results_dir}')

    cell_state_results_dir = analysis_dirs['cell_state_prediction']

    manifest_filename = os.path.join(cell_state_results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

    filenames = scgenome.loaders.utils.find_filenames(manifest['filenames'], 'cell_state_prediction.csv')

    if len(filenames) != 1:
        raise ValueError(f'found {len(filenames)} filenames with suffix cell_state_prediction.csv')

    filepath = os.path.join(cell_state_results_dir, filenames[0])
    data = pd.read_csv(filepath)

    return data


def load_qc_data(
        results_dir,
        sample_ids=None,
        subsample=None,
        additional_hmmcopy_reads_cols=None,
    ):
    """ Load qc data (align, hmmcopy, annotation)
    
    Args:
        results_dir (str): results directory to load from.
        sample_ids (list of str, optional): Set of sample ids to filter for. Defaults to None.
        subsample (float, optional): Proportion of the cells to downsample to. Defaults to None.
        additional_hmmcopy_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    """

    ticket_results_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    results_tables = load_align_data(ticket_results_dirs['align'])
    hmmcopy_results_tables = load_hmmcopy_data(
        ticket_results_dirs['hmmcopy'],
        additional_reads_cols=additional_hmmcopy_reads_cols)

    results_tables.update(hmmcopy_results_tables)

    # Load annotation tables if they exist otherwise create merge of hmmcopy/align
    if 'annotation' in ticket_results_dirs:
        annotation_results_tables = load_annotation_data(ticket_results_dirs['annotation'])
        results_tables['annotation_metrics'] = annotation_results_tables['annotation_metrics']

    else:
        results_tables['annotation_metrics'] = _calculate_annotation_metrics(results_tables)

    # For older results annotation metrics will not contain s phase, load directly
    if 'is_s_phase' not in results_tables['annotation_metrics']:
        cell_state = load_cell_state_prediction(ticket_results_dirs['cell_state_prediction'])
        results_tables['annotation_metrics'] = results_tables['annotation_metrics'].merge(
            cell_state[['cell_id', 'is_s_phase', 'is_s_phase_prob']].drop_duplicates())

    # Optionally select specific samples
    if sample_ids is not None:
        for table_name, table_data in results_tables.items():
            results_tables[table_name] = table_data[table_data['sample_id'].isin(sample_ids)]

    # Optionally subsample cells
    if subsample is not None:
        cell_ids = (
            results_tables['hmmcopy_metrics'][['cell_id']]
            .drop_duplicates().sample(frac=subsample))
        for table_name in results_tables.keys():
            results_tables[table_name] = results_tables[table_name].merge(cell_ids)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables

