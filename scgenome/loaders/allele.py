import logging
import yaml
import os
import pandas as pd

import scgenome.utils
import scgenome.loaders.utils
import scgenome.csvutils


def load_haplotype_allele_data(
        results_dir,
    ):
    """ Load the haplotype allele count data from the pseudobulk results paths
    
    Args:
        pseudobulk_dir (str): results directory

    Returns:
        dict of pandas.DataFrame: Haplotype allele data
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['pseudobulk']
        suffix = 'allele_counts.csv'

    elif 'count_haps' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['count_haps']
        suffix = 'allele_counts.tsv'

    else:
        raise ValueError(f'no pseudobulk found for directory {results_dir}')

    allele_counts = []

    files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, suffix)

    for sample_id, library_id, filepath in files:
        logging.info('Loading haplotype allele counts from {}'.format(filepath))

        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv(
            dtypes_override={
                'chromosome': 'category',
                'cell_id': 'category',
        })

        if library_id is not None:
            data['library_id'] = pd.Series([library_id], dtype="category")

        if sample_id is not None:
            data['sample_id'] = pd.Series([sample_id], dtype="category")

        logging.info(f'Loaded haplotype allele counts table with shape {data.shape}, memory {data.memory_usage().sum()}')

        allele_counts.append(data)

    allele_counts = scgenome.utils.concat_with_categories(allele_counts, ignore_index=True)

    logging.info(f'Loaded all haplotype allele counts table with shape {allele_counts.shape}, memory {allele_counts.memory_usage().sum()}')

    return {
        'allele_counts': allele_counts,
    }

