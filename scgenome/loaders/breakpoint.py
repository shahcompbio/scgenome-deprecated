import pandas as pd
import wgs_analysis.algorithms.rearrangement
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils


def load_destruct_from_results(results_dir):
    """ Load destruct breakpoint data from a results directory.

    Args:
        results_dir (str): results directory
    """

    destruct_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    return load_destruct(destruct_filepath)


def load_destruct(destruct_filepath):
    """ Read destruct breakpoint data from a file.

    Args:
        destruct_filepath (str): results filepath
    """

    csv_input = scgenome.csvutils.CsvInput(destruct_filepath)
    data = csv_input.read_csv()

    data['chromosome_1'] = data['chromosome_1'].astype(str)
    data['chromosome_2'] = data['chromosome_2'].astype(str)

    return data


def load_destruct_counts_from_results(results_dir):
    """ Load destruct breakpoint count data from a results directory.

    Args:
        results_dir (str): results directory
    """

    destruct_counts_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'destruct_cell_counts.csv.gz', analysis_type='breakpoint_calling')

    return load_destruct_counts(destruct_counts_filepath)


def load_destruct_counts(destruct_counts_filepath):
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


def load_lumpy_from_results(results_dir, standardize_columns=False):
    """ Load lumpy breakpoint data from a results directory.

    Args:
        results_dir(str): results directory
        
    KwArgs:
        standardize_columns(bool): rename columns to be similar to destruct
    """

    lumpy_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'lumpy_breakpoints.csv.gz', analysis_type='breakpoint_calling')

    return load_lumpy(lumpy_filepath, standardize_columns=standardize_columns)


def load_lumpy(lumpy_filepath, standardize_columns=False):
    """ Load lumpy breakpoint data from a file.

    Args:
        lumpy_filepath(str): results filepath
        
    KwArgs:
        standardize_columns(bool): rename columns to be similar to destruct
    """

    csv_input = scgenome.csvutils.CsvInput(lumpy_filepath)
    data = csv_input.read_csv()

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


def load_lumpy_counts_from_results(results_dir):
    """ Load lumpy breakpoint counts data from a results directory.

    Args:
        results_dir (str): results directory
    """

    lumpy_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'lumpy_breakpoints_evidence.csv.gz', analysis_type='breakpoint_calling')

    return load_lumpy_counts(lumpy_filepath)


def load_lumpy_counts(lumpy_filepath):
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


def load_consensus_from_results(results_dir):
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

    destruct_breakpoints = scgenome.loaders.breakpoint.load_destruct(destruct_filepath)
    lumpy_breakpoints = scgenome.loaders.breakpoint.load_lumpy(lumpy_filepath, standardize_columns=True)

    destruct_lumpy_matches = wgs_analysis.algorithms.rearrangement.match_breakpoints(
        destruct_breakpoints, lumpy_breakpoints)

    # Filter destruct by lumpy matches
    destruct_breakpoints = destruct_breakpoints.merge(
        destruct_lumpy_matches[['reference_prediction_id']].drop_duplicates()
        .rename(columns={'reference_prediction_id': 'prediction_id'}))

    destruct_breakpoint_counts = scgenome.loaders.breakpoint.load_destruct_counts(destruct_counts_filepath)

    return {
        'breakpoint_data': destruct_breakpoints,
        'breakpoint_count_data': destruct_breakpoint_counts,
    }


# Deprecated from here
# VVV

def load_breakpoint_annotation_data(files, is_lumpy=False):
    """ Load breakpoint data from a pseudobulk run.

    Args:
        files (str): results directory
        suffix (str): suffix of breakpoint annotation tables
    """

    chrom_1_colname = "chromosome_1"
    chrom_2_colname = "chromosome_2"

    if is_lumpy:
        chrom_1_colname = "chrom1"
        chrom_2_colname = "chrom2"

    breakpoint_data = []

    for sample_id, library_id, filepath in files:
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()

        data[chrom_1_colname] = data[chrom_1_colname].astype(str)
        data[chrom_2_colname] = data[chrom_2_colname].astype(str)

        if library_id is not None:
            data['library_id'] = library_id

        if sample_id is not None:
            data['sample_id'] = sample_id

        breakpoint_data.append(data)

    if len(breakpoint_data) == 0:
        return pd.DataFrame(), pd.DataFrame()

    breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

    return breakpoint_data


def load_breakpoint_count_data(files):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        pseudobulk_dir (str): results directory
        suffix (str): suffix of breakpoint count tables
    """
    breakpoint_count_data = []

    for sample_id, library_id, filepath in files:
        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv()

        if library_id is not None:
            data['library_id'] = pd.Series([library_id], dtype="category")

        if sample_id is not None:
            data['sample_id'] = pd.Series([sample_id], dtype="category")

        breakpoint_count_data.append(data)

    breakpoint_count_data = pd.concat(breakpoint_count_data, ignore_index=True)
    breakpoint_count_data = breakpoint_count_data.rename(columns={'cluster_id': 'prediction_id'})

    # KLUDGE: normal reads are not filtered properly, filter by their prefix, and having '-' in cell id
    breakpoint_count_data = breakpoint_count_data.loc[~breakpoint_count_data['cell_id'].str.startswith('HS'), :]
    breakpoint_count_data = breakpoint_count_data.loc[breakpoint_count_data['cell_id'].apply(lambda a: '-' in a), :]

    return breakpoint_count_data


def load_breakpoint_data_from_files(annotation_file, counts_file, lumpy=False):
    annotation_files = scgenome.loaders.utils._prep_filenames_for_loading(annotation_file)

    breakpoint_data = load_breakpoint_annotation_data(annotation_files, is_lumpy=lumpy)

    count_files = scgenome.loaders.utils._prep_filenames_for_loading(counts_file)

    breakpoint_count_data = load_breakpoint_count_data(count_files)

    return _process_breakpoint_data(breakpoint_data, breakpoint_count_data)


def load_breakpoint_data(
        results_dir,
):
    """ Load breakpoint count data from a pseudobulk run.

    Args:
        results_dir (str): results directory to load from.
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        breakpoint_calling_dir = analysis_dirs['pseudobulk']
        annotation_suffix = 'destruct.csv.gz'
        count_suffix = 'cell_counts_destruct.csv.gz'

    elif 'breakpoint_calling' in analysis_dirs:
        breakpoint_calling_dir = analysis_dirs['breakpoint_calling']

        if len(breakpoint_calling_dir) == 0:
            raise ValueError(f'found {len(breakpoint_calling_dir)} dirs for breakpoint_calling')

        elif len(breakpoint_calling_dir) > 1:
            if filter_sample_id is None:
                raise ValueError(f'found {len(breakpoint_calling_dir)} without filter_sample_id')

            filtered_breakpoint_calling_dir = list(
                filter(lambda a: f'sample_{filter_sample_id}' in a, breakpoint_calling_dir))

            if len(filtered_breakpoint_calling_dir) != 1:
                raise ValueError(
                    f'found {len(filtered_breakpoint_calling_dir)} in {breakpoint_calling_dir} matching filter_sample_id')

            breakpoint_calling_dir = filtered_breakpoint_calling_dir

        breakpoint_calling_dir = breakpoint_calling_dir[0]
        annotation_suffix = 'destruct_breakpoints.csv.gz'
        count_suffix = 'destruct_cell_counts.csv.gz'

    else:
        raise ValueError(f'no breakpoints found for directory {results_dir}')

    annotation_files = scgenome.loaders.utils.get_pseudobulk_files(
        breakpoint_calling_dir, annotation_suffix)

    breakpoint_data = load_breakpoint_annotation_data(annotation_files)

    count_files = scgenome.loaders.utils.get_pseudobulk_files(
        breakpoint_calling_dir, count_suffix)

    breakpoint_count_data = load_breakpoint_count_data(count_files)

    return _process_breakpoint_data(breakpoint_data, breakpoint_count_data)


def _process_breakpoint_data(breakpoint_data, breakpoint_count_data):
    # TODO: fix upstream
    for col in ('prediction_id', 'position_1', 'position_2', 'read_count'):
        for df in (breakpoint_data, breakpoint_count_data):
            if col in df:
                df[col] = df[col].astype(int)

    return {
        'breakpoint_data': breakpoint_data,
        'breakpoint_count_data': breakpoint_count_data,
    }
