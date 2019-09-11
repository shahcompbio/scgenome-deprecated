import logging
import yaml
import os
import pandas as pd

import scgenome.utils
import scgenome.loaders.utils


def load_haplotype_allele_data(results_dir):
    """ Load the haplotype allele count data from the pseudobulk results paths
    
    Args:
        results_dir (str): results directory to load from.
    
    Returns:
        dict of pandas.DataFrame: Haplotype allele data
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' not in analysis_dirs:
        raise ValueError(f'no pseudobulk found for directory {results_dir}')

    pseudobulk_dir = analysis_dirs['pseudobulk']

    allele_counts = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, 'allele_counts.csv'):
        logging.info('Loading haplotype allele counts from {}'.format(filepath))

        # HACK: temp fix until the pipeline outputs csv headers
        names = None
        firstline = open(filepath).readline()
        if 'chromosome' not in firstline:
            names = [
                'start',
                'end',
                'hap_label',
                'allele_id',
                'readcount',
                'chromosome',
                'cell_id',
            ]

        data = pd.read_csv(
            filepath,
            names=names,
            dtype={
                'chromosome': 'category',
                'cell_id': 'category',
            })

        logging.info('Loaded haplotype allele counts table with shape {}'.format(data.shape))

        allele_counts.append(data)

    allele_counts = scgenome.utils.concat_with_categories(allele_counts, ignore_index=True)

    logging.info('Loaded all haplotype allele counts table with shape {}'.format(allele_counts.shape))

    return {
        'allele_counts': allele_counts,
    }

