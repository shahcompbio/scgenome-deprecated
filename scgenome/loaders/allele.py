import logging

import pandas as pd
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils


def load_haplotype_allele_counts_results(results_dir):
    """ Load haplotype allele counts from a results directory.

    Args:
        results_dir (str): results directory
    """

    allele_counts_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'allele_counts.csv.gz', analysis_type='count_haps')

    return load_haplotype_allele_counts_files(allele_counts_filepath)


def load_haplotype_allele_counts_files(allele_counts_filepath):
    """ Load haplotype allele counts from a file.

    Args:
        allele_counts_filepath (str): results filepath
    """

    csv_input = scgenome.csvutils.CsvInput(allele_counts_filepath)
    data = csv_input.read_csv(
        dtypes_override={
            'chromosome': 'category',
            'cell_id': 'category',
    })

    return data

