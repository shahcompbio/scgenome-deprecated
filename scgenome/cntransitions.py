from numba import jit
import numpy as np
import pandas as pd


@jit(nopython=True)
def count_pairs(group_items):
    num_items = group_items.shape[1]
    pair_counts = np.zeros((num_items, num_items))
    for group_idx in range(group_items.shape[0]):
        items = np.where(group_items[group_idx] == 1)[0]
        for item_idx1 in items:
            for item_idx2 in items:
                pair_counts[item_idx1, item_idx2] += 1
    return pair_counts


def generate_cn_transitions(cn_data):
    """ Generate a list of copy number transition regions per cell.

    Args:
        cn_data (DataFrame): copy number data
    
    Returns:
        DataFrame: copy number transitions per cell

    """

    # Build a list of transition bins
    cn_transitions = cn_data.sort_values(['cell_id', 'chr', 'start'])[['cell_id', 'chr', 'start', 'end', 'state']]
    cn_transitions['state_diff'] = cn_transitions['state'].diff().fillna(0).astype(int)
    cn_transitions['chr_diff'] = (cn_transitions['chr'].shift(1) == cn_transitions['chr']) * 1
    cn_transitions['transition'] = (cn_transitions['state_diff'] != 0) & (cn_transitions['chr_diff'] == 1)
    cn_transitions['start'] -= int(500000/2)
    cn_transitions['end'] -= int(500000/2)
    cn_transitions['orientation'] = '+'
    cn_transitions.loc[cn_transitions['state_diff'] > 0, 'orientation'] = '-'
    cn_transitions = cn_transitions.query('transition')
    cn_transitions = cn_transitions.merge(
        cn_transitions[['chr', 'start', 'end', 'orientation']]
            .drop_duplicates()
            .reset_index(drop=True)
            .reset_index()
            .rename(columns={'index': 'region_index'}))

    return cn_transitions


def generate_cn_transition_counts(cn_data):
    """ Generate a list of copy number transition regions and their cell counts.

    Args:
        cn_data (DataFrame): copy number data
    
    Returns:
        DataFrame: copy number transitions cell counts

    """

    cn_transitions = generate_cn_transitions(cn_data)

    # Build a table of transition cell counts
    cn_transition_counts = (
        cn_transitions
            .groupby(['chr', 'start', 'end', 'region_index', 'orientation'], observed=True)
            .size().rename('cell_count').reset_index())

    return cn_transition_counts


def generate_cn_transition_pair_counts(cn_data):
    """ Generate a list of copy number transition pairs and their cell counts.

    Args:
        cn_data (DataFrame): copy number data
    
    Returns:
        DataFrame: copy number transitions pairs and counts of supporting cells

    """

    cn_transitions = generate_cn_transitions(cn_data)

    # Build a matrix of pairs of transitions in the same cell
    cell_transitions_matrix = pd.crosstab(
        index=cn_transitions['cell_id'],
        columns=cn_transitions['region_index'],
    )

    # Possibly not necessary, but ensure the column ordering matches the column indexing
    cell_transitions_matrix = cell_transitions_matrix.reindex(
        columns=range(max(cell_transitions_matrix.columns) + 1))

    # For each pair of transitions, count the number of cells with both
    # returned as a transition by transition matrix
    cell_region_counts = count_pairs(cell_transitions_matrix.values)

    # Create a stacked format matrix that does not include 0 entries
    cell_region_counts[cell_region_counts == 0] = np.NaN

    cell_region_counts = pd.DataFrame(
        cell_region_counts,
        columns=cell_transitions_matrix.columns,
        index=cell_transitions_matrix.columns,
    )

    cn_transitions_matrix = (
        cell_region_counts
            .stack().rename('pair_cell_count')
            .rename_axis(index=['region_index_1', 'region_index_2'])
            .reset_index())

    cn_transitions_matrix['pair_cell_count'] = cn_transitions_matrix['pair_cell_count'].astype(int)

    return cn_transitions_matrix
