import os
import yaml
import pandas as pd

import scgenome.loaders.utils
from scgenome.loaders.align import load_align_data
from scgenome.loaders.align import load_align_data_from_files

from scgenome.loaders.hmmcopy import load_hmmcopy_data
from scgenome.loaders.hmmcopy import load_hmmcopy_data_from_filename

from scgenome.loaders.annotation import load_annotation_data
from scgenome.loaders.annotation import load_annotation_data_from_file


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

    if len(cell_state_results_dir) != 1:
        raise ValueError(f"found {len(cell_state_results_dir)} directories with cell_state_prediction results")
    cell_state_results_dir = cell_state_results_dir[0]

    manifest_filename = os.path.join(cell_state_results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

    filenames = scgenome.loaders.utils.find_filenames(manifest['filenames'], 'cell_state_prediction.csv')

    if len(filenames) != 1:
        raise ValueError(f'found {len(filenames)} filenames with suffix cell_state_prediction.csv')

    filepath = os.path.join(cell_state_results_dir, filenames[0])
    data = pd.read_csv(filepath)

    if 'is_s_phase' in data:
        data['is_s_phase'] = data['is_s_phase'].fillna(False).astype(bool)

    return data


def load_qc_data_from_files(hmmcopy_reads, hmmcopy_segs, 
    hmmcopy_metrics, alignment_metrics, gc_metrics, annotation_metrics=None, 
    sample_id=None, additional_hmmcopy_reads_cols=None,
):

    results_tables = load_align_data_from_files(alignment_metrics, 
        gc_metrics=gc_metrics
    )

    hmmcopy_results_tables = load_hmmcopy_data_from_filename(hmmcopy_reads, 
        hmmcopy_segs, hmmcopy_metrics,
        additional_reads_cols=additional_hmmcopy_reads_cols
    )
    results_tables.update(hmmcopy_results_tables)

    if annotation_metrics:
        annotation_results_tables = load_annotation_data_from_file(annotation_metrics)
    else:
        annotation_results_tables = _calculate_annotation_metrics(results_tables)
    results_tables.update(annotation_results_tables)

  

    if sample_id is not None:
        results_tables = _sample_id_filter(results_tables, sample_id)

    
    scgenome.utils.union_categories(results_tables.values())
    return results_tables

def load_qc_data(
        results_dir,
        sample_ids=None,
        additional_hmmcopy_reads_cols=None,
    ):
    """ Load qc data (align, hmmcopy, annotation)
    
    Args:
        results_dir (str): results directory to load from.
        sample_ids (list of str, optional): Set of sample ids to filter for. Defaults to None.
        additional_hmmcopy_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    """
    
    ticket_results_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if len(ticket_results_dirs['align']) != 1:
        raise ValueError(f"found {len(ticket_results_dirs['align'])} directories with align results")

    results_tables = load_align_data(ticket_results_dirs['align'][0])

    if len(ticket_results_dirs['hmmcopy']) != 1:
        raise ValueError(f"found {len(ticket_results_dirs['hmmcopy'])} directories with hmmcopy results")

    hmmcopy_results_tables = load_hmmcopy_data(
        ticket_results_dirs['hmmcopy'][0],
        additional_reads_cols=additional_hmmcopy_reads_cols)

    results_tables.update(hmmcopy_results_tables)

    # Load annotation tables if they exist otherwise create merge of hmmcopy/align
    if 'annotation' in ticket_results_dirs:
        if len(ticket_results_dirs['annotation']) != 1:
            raise ValueError(f"found {len(ticket_results_dirs['annotation'])} directories with annotation results")

        annotation_results_tables = load_annotation_data(ticket_results_dirs['annotation'][0])

        results_tables['annotation_metrics'] = annotation_results_tables['annotation_metrics']

    else:
        results_tables['annotation_metrics'] = _calculate_annotation_metrics(results_tables)

    # For older results annotation metrics will not contain s phase, load directly
    if 'is_s_phase' not in results_tables['annotation_metrics']:
        if 'cell_state_prediction' not in ticket_results_dirs:
            raise ValueError(f'no cell state predictions found in {results_dir}')
        if len(ticket_results_dirs['cell_state_prediction']) != 1:
            raise ValueError(f"found {len(ticket_results_dirs['cell_state_prediction'])} directories with cell_state_prediction results")
        cell_state = load_cell_state_prediction(ticket_results_dirs['cell_state_prediction'][0])
        results_tables['annotation_metrics'] = results_tables['annotation_metrics'].merge(
            cell_state[['cell_id', 'is_s_phase', 'is_s_phase_prob']].drop_duplicates())

    if sample_ids is not None:
        results_tables = _sample_id_filter(results_tables, sample_ids)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def _sample_id_filter(results_tables, sample_ids):
    # Optionally select specific samples
    for table_name, table_data in results_tables.items():
        results_tables[table_name] = table_data[table_data['sample_id'].isin(sample_ids)]

    return results_tables
