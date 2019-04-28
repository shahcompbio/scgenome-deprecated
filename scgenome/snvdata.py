import os
import logging
import pandas as pd
import matplotlib.pyplot as plt

from datamanagement.miscellaneous.hdf5helper import read_python2_hdf5_dataframe
import scgenome.utils


def get_highest_snpeff_effect(snpeff_data):
    """ Select the highest ranked effect from snpeff data.
    """
    ordered_effect_impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    snpeff_data = snpeff_data.merge(
        pd.DataFrame({
            'effect_impact': ordered_effect_impacts,
            'effect_impact_rank': range(len(ordered_effect_impacts))}))
    snpeff_data = (
        snpeff_data.sort_values(['chrom', 'coord', 'ref', 'alt', 'effect_impact_rank'], ascending=True)
        .groupby(['chrom', 'coord', 'ref', 'alt'], sort=False)
        .first()
        .reset_index())
    snpeff_data = snpeff_data[[
        'chrom', 'coord', 'ref', 'alt',
        'gene_name', 'effect', 'effect_impact', 'amino_acid_change',
    ]]
    return snpeff_data


categorical_columns = [
    'chrom',
    'ref',
    'alt',
    'gene_name',
    'effect',
    'effect_impact',
    'amino_acid_change',
    'tri_nucleotide_context',
]

default_museq_filter = 0.9
default_strelka_filter = 20.

def get_snv_results(dest, museq_filter=None, strelka_filter=None):
    """ Collate snv results into a single table.
    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    logging.info('starting load')
    mappability = read_python2_hdf5_dataframe(dest, '/snv/mappability')
    mappability = mappability[mappability['mappability'] > 0.99]
    mappability['chrom'] = mappability['chrom'].astype(str)

    strelka_results = read_python2_hdf5_dataframe(dest, '/strelka/vcf').rename(columns={'score': 'strelka_score'})
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = read_python2_hdf5_dataframe(dest, '/museq/vcf').rename(columns={'score': 'museq_score'})
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = read_python2_hdf5_dataframe(dest,'/snv/cosmic_status')
    logging.info('cosmic table with shape {}'.format(cosmic.shape))
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = read_python2_hdf5_dataframe(dest,'/snv/snpeff')
    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info('snpeff table with shape {}'.format(snpeff.shape))

    tnc = read_python2_hdf5_dataframe(dest,'/snv/tri_nucleotide_context')

    # data = read_python2_hdf5_dataframe(dest, '/snv_allele_counts')
    # logging.info('total snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = mappability#data.merge(mappability)
    logging.info('post mappability with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(strelka_results, how='left')
    logging.info('post strelka with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(museq_results, how='left')
    logging.info('post museq with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(cosmic, how='left')
    logging.info('post cosmic with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(snpeff, how='left')
    logging.info('post snpeff with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data.merge(tnc, how='left')

    data = data[data['museq_score'] > museq_filter]
    logging.info('post museq filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    data = data[data['strelka_score'] > strelka_filter]
    logging.info('post strelka filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    logging.info('finishing load with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))

    for column in categorical_columns:
        data[column] = data[column].astype('category')

    return data


def load_snv_annotation_data(dataset_filepaths, museq_filter=None, strelka_filter=None):
    snv_data = []
    for filepath in dataset_filepaths:
        if filepath.endswith('_snv_annotations.h5'):
            snv_data.append(get_snv_results(
                filepath,
                museq_filter=museq_filter,
                strelka_filter=strelka_filter))

    snv_data = scgenome.utils.concat_with_categories(snv_data, ignore_index=True)

    return snv_data


def load_snv_count_data(dataset_filepaths, positions):
    snv_count_data = []

    for filepath in dataset_filepaths:
        if filepath.endswith('_snv_union_counts.csv.gz'):
            logging.info('Loading snv counts from {}'.format(filepath))

            data = pd.read_csv(
                filepath,
                dtype={
                    'chrom': 'category',
                    'ref': 'category',
                    'alt': 'category',
                    'cell_id': 'category',
                })

            logging.info('Loaded snv counts table with shape {}'.format(data.shape))

            scgenome.utils.union_categories(
                data, positions,
                cols=['chrom', 'ref', 'alt'])
            data = data.merge(positions, how='inner')

            logging.info('Filtered snv counts table to shape {}'.format(data.shape))

            snv_count_data.append(data)

    snv_count_data = scgenome.utils.concat_with_categories(snv_count_data, ignore_index=True)

    logging.info('Loaded all snv counts tables with shape {}'.format(snv_count_data.shape))

    return snv_count_data


def load_snv_data(
        dataset_filepaths,
        museq_filter=None,
        strelka_filter=None,
        num_cells_threshold=None,
        sum_alt_threshold=None,
        figures_prefix=None,
):
    snv_data = load_snv_annotation_data(
        dataset_filepaths,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter)

    assert not snv_data['coord'].isnull().any()

    positions = snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

    snv_count_data = load_snv_count_data(dataset_filepaths, positions)
    snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
    snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

    assert not snv_count_data['alt_counts'].isnull().any()

    # Calculate cell counts
    cell_counts = (
        snv_count_data
        .query('alt_counts > 0')
        .drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
        .groupby(['chrom', 'coord', 'ref', 'alt']).size().rename('num_cells')
        .reset_index())

    if figures_prefix is not None:
        fig = plt.figure(figsize=(4, 4))
        cell_counts['num_cells'].hist(bins=50)
        fig.savefig(figures_prefix + '_snv_cell_counts.pdf', bbox_inches='tight')

    snv_data = snv_data.merge(cell_counts, how='left')
    assert not snv_data['num_cells'].isnull().any()

    # Calculate total alt counts for each SNV
    sum_alt_counts = (
        snv_count_data
        .groupby(['chrom', 'coord', 'ref', 'alt'])['alt_counts']
        .sum().rename('sum_alt_counts')
        .reset_index())

    if figures_prefix is not None:
        fig = plt.figure(figsize=(4, 4))
        sum_alt_counts['sum_alt_counts'].hist(bins=50)
        fig.savefig(figures_prefix + '_snv_alt_counts.pdf', bbox_inches='tight')

    # Filter SNVs by num cells in which they are detected
    if num_cells_threshold is not None:
        filtered_cell_counts = cell_counts.query(
            'num_cells >= {}'.format(num_cells_threshold))
        logging.info('Filtering {} of {} SNVs by num_cells >= {}'.format(
            len(filtered_cell_counts.index), len(cell_counts.index), num_cells_threshold))
        snv_data = snv_data.merge(filtered_cell_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    # Filter SNVs by total alt counts
    if sum_alt_threshold is not None:
        filtered_sum_alt_counts = sum_alt_counts.query(
            'sum_alt_counts >= {}'.format(sum_alt_threshold))
        logging.info('Filtering {} of {} SNVs by sum_alt_counts >= {}'.format(
            len(filtered_sum_alt_counts.index), len(sum_alt_counts.index), sum_alt_threshold))
        snv_data = snv_data.merge(filtered_sum_alt_counts[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    snv_count_data = snv_count_data.merge(snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates())

    return snv_data, snv_count_data



