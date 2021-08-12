import pandas as pd
import wgs_analysis.algorithms.rearrangement
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils


def load_destruct_results(results_dir):
    """ Load destruct breakpoint data from a results directory.

    Args:
        results_dir (str): results directory
    """

    destruct_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    return load_destruct_files(destruct_filepath)


def load_destruct_files(destruct_filepath):
    """ Read destruct breakpoint data from a file.

    Args:
        destruct_filepath (str): results filepath
    """

    csv_input = scgenome.csvutils.CsvInput(destruct_filepath)
    data = csv_input.read_csv()

    data['chromosome_1'] = data['chromosome_1'].astype(str)
    data['chromosome_2'] = data['chromosome_2'].astype(str)

    return data


def load_destruct_counts_results(results_dir):
    """ Load destruct breakpoint count data from a results directory.

    Args:
        results_dir (str): results directory
    """

    destruct_counts_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_cell_counts.csv.gz', analysis_type='breakpoint_calling')

    return load_destruct_counts_files(destruct_counts_filepath)


def load_destruct_counts_files(destruct_counts_filepath):
    """ Load destruct breakpoint count data from a file.

    Args:
        destruct_counts_filepath (str): results filepath
    """

    csv_input = scgenome.csvutils.CsvInput(destruct_counts_filepath)
    data = csv_input.read_csv()

    data = data.rename(columns={'cluster_id': 'prediction_id'})

    # KLUDGE: normal reads are not filtered properly, filter by their prefix, and having '-' in cell id
    data = data.loc[~data['cell_id'].str.startswith('HS'), :]
    data = data.loc[data['cell_id'].apply(lambda a: '-' in a), :]
    
    return data


def load_lumpy_results(results_dir, standardize_columns=False):
    """ Load lumpy breakpoint data from a results directory.

    Args:
        results_dir(str): results directory
        
    KwArgs:
        standardize_columns(bool): rename columns to be similar to destruct
    """

    lumpy_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'lumpy_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    return load_lumpy_files(lumpy_filepath, standardize_columns=standardize_columns)


def load_lumpy_files(lumpy_filepath, standardize_columns=False):
    """ Load lumpy breakpoint data from a file.

    Args:
        lumpy_filepath(str): results filepath

    KwArgs:
        standardize_columns(bool): rename columns to be similar to destruct
    """

    csv_input = scgenome.csvutils.CsvInput(lumpy_filepath)
    data = csv_input.read_csv()

    count_cols = ['normal_PE', 'tumour_SR', 'normal_SR', 'tumour_PE']
    data[count_cols] = data[count_cols].fillna(0).astype(int)

    data = data[(data['normal_PE'] == 0) & (data['normal_SR'] == 0)]
    
    data['chrom1'] = data['chrom1'].astype(str)
    data['chrom2'] = data['chrom2'].astype(str)

    # Standardize columns
    if standardize_columns:
        data = data.rename(columns={
            'breakpoint_id': 'prediction_id',
            'chrom1': 'chromosome_1',
            'chrom2': 'chromosome_2',
            'strand1': 'strand_1',
            'strand2': 'strand_2',
        })

        type_translate = {
            'DELETION': 'deletion',
            'DUPLICATION': 'duplication',
            'INVERSION': 'inversion',
            'INTERCHROM': 'translocation',
        }
        data['type'] = data['type'].apply(lambda a: type_translate[a])

        # Note: its unclear whether lumpy uses 0-based or 1-based coordinates
        data['position_1'] = data[['confidence_interval_start1', 'confidence_interval_end1']].mean(axis=1)
        data['position_2'] = data[['confidence_interval_start2', 'confidence_interval_end2']].mean(axis=1)

        data['position_1'] = data['position_1'].astype(int)
        data['position_2'] = data['position_2'].astype(int)

    return data


def load_lumpy_counts_results(results_dir):
    """ Load lumpy breakpoint counts data from a results directory.

    Args:
        results_dir (str): results directory
    """

    lumpy_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'lumpy_breakpoints_evidence.csv.gz', analysis_type='breakpoint_calling')

    return load_lumpy_counts_files(lumpy_filepath)


def load_lumpy_counts_files(lumpy_filepath):
    """ Load lumpy breakpoint counts data from a file.

    Args:
        lumpy_filepath (str): results filepath
    """

    csv_input = scgenome.csvutils.CsvInput(lumpy_filepath)
    data = csv_input.read_csv()

    # KLUDGE: normal reads are not filtered properly, filter by their prefix, and having '-' in cell id
    data = data.loc[~data['cell_id'].str.startswith('HS'), :]
    data = data.loc[data['cell_id'].apply(lambda a: '-' in a), :]

    return data


def load_consensus_breakpoint_results(results_dir):
    """ Load / calculate consensus breakpoints from a results directory.

    Args:
        results_dir (str): results directory
    """

    destruct_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    lumpy_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'lumpy_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    destruct_counts_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_cell_counts.csv.gz', analysis_type='breakpoint_calling')

    return load_consensus_breakpoint_files(destruct_filepath, lumpy_filepath, destruct_counts_filepath)


def load_consensus_breakpoint_files(destruct_filepath, lumpy_filepath, destruct_counts_filepath):
    """ Load / calculate consensus breakpoints from files.

    Args:
        results_dir (str): results directory
    """

    destruct_breakpoints = scgenome.loaders.breakpoint.load_destruct_files(destruct_filepath)
    lumpy_breakpoints = scgenome.loaders.breakpoint.load_lumpy_files(lumpy_filepath, standardize_columns=True)

    destruct_breakpoint_counts = scgenome.loaders.breakpoint.load_destruct_counts_files(destruct_counts_filepath)

    destruct_lumpy_matches = wgs_analysis.algorithms.rearrangement.match_breakpoints(
        lumpy_breakpoints, destruct_breakpoints, window_size=200)

    matched_destruct_predictions = (
        destruct_lumpy_matches[['target_id']].drop_duplicates()
        .rename(columns={'target_id': 'prediction_id'}))

    # Filter destruct by lumpy matches
    destruct_breakpoints = destruct_breakpoints.merge(matched_destruct_predictions)
    destruct_breakpoint_counts = destruct_breakpoint_counts.merge(matched_destruct_predictions)

    return {
        'breakpoint_data': destruct_breakpoints,
        'breakpoint_count_data': destruct_breakpoint_counts,
    }

