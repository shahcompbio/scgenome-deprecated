import logging
import yaml
import os
import pandas as pd
import numpy as np

import scgenome.utils
import scgenome.loaders.utils
import scgenome.csvutils


default_museq_filter = 0.9
default_strelka_filter = 20.


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


def load_snv_count_data(pseudobulk_dir, suffix, positions, filter_sample_id=None, filter_library_id=None):
    """ Load per cell SNV count data
    
    Args:
        pseudobulk_dir (str): pseudobulk results directory
        suffix (str): suffix of snv count tables
        positions (pandas.DataFrame): restrict to the specified positions

    Kwargs:
        filter_sample_id (str): restrict to specific sample id
        filter_library_id (str): restrict to specific library id
    
    Returns:
        pandas.DataFrame: SNV alt and ref counts per cell
    """
    snv_count_data = []

    files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, suffix)

    for sample_id, library_id, filepath in files:
        logging.info('Loading snv counts from {}'.format(filepath))

        if sample_id is not None and filter_sample_id is not None and sample_id != filter_sample_id:
            logging.info('skipping {sample_id}', sample_id, filter_sample_id)
            continue

        if library_id is not None and filter_sample_id is not None and library_id != filter_library_id:
            logging.info('skipping', library_id, filter_library_id)
            continue

        data = []
        csv_input = scgenome.csvutils.CsvInput(filepath)

        chunk_iter = csv_input.read_csv(
            chunksize=10**6,
            dtypes_override={
                'chrom': 'category',
                'ref': 'category',
                'alt': 'category',
                'cell_id': 'category',
                'sample_id': 'category',
                'library_id': 'category',
            },
        )

        for chunk in chunk_iter:
            scgenome.utils.union_categories(
                [chunk, positions],
                cat_cols=['chrom', 'ref', 'alt'])
            chunk = chunk.merge(positions)

            if filter_sample_id is not None and 'sample_id' in chunk:
                chunk = chunk[chunk['sample_id'] == filter_sample_id]

            if filter_library_id is not None and 'library_id' in chunk:
                chunk = chunk[chunk['library_id'] == filter_library_id]

            data.append(chunk)

        data = scgenome.utils.concat_with_categories(data, ignore_index=True)

        if library_id is not None:
            data['library_id'] = pd.Series([library_id], dtype="category")

        if sample_id is not None:
            data['sample_id'] = pd.Series([sample_id], dtype="category")

        logging.info(f'Loaded snv counts table with shape {data.shape}, memory {data.memory_usage().sum()}')

        scgenome.utils.union_categories(
            [data, positions],
            cat_cols=['chrom', 'ref', 'alt'])
        data = data.merge(positions, how='inner')

        logging.info(f'Filtered snv counts table to shape {data.shape}, memory {data.memory_usage().sum()}')

        snv_count_data.append(data)

    snv_count_data = scgenome.utils.concat_with_categories(snv_count_data, ignore_index=True)

    logging.info(f'Loaded all snv counts tables with shape {snv_count_data.shape}, memory {snv_count_data.memory_usage().sum()}')

    return snv_count_data



def load_snv_annotation_table(pseudobulk_dir, table_name):
    """ Load SNV annotation data

    Args:
        pseudobulk_dir (str): pseudobulk results directory
        table_name (str): name of annotation table to load.

    Returns:
        pandas.DataFrame: SNVs annotation data per sample / library
    """
    snv_data = []

    files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, f'snv_{table_name}.csv.gz')

    for sample_id, library_id, filepath in files:
        logging.info(f'Loading snv {table_name} annotations from {filepath}')

        csv_input = scgenome.csvutils.CsvInput(filepath)
        data = csv_input.read_csv(
            dtypes_override={
                'chrom': 'category',
                'ref': 'category',
                'alt': 'category',
                'cell_id': 'category',
                'effect': 'category',
                'effect_impact': 'category',
                'functional_class': 'category',
                'codon_change': 'category',
                'amino_acid_change': 'category',
                'gene_name': 'category',
                'transcript_biotype': 'category',
                'gene_coding': 'category',
                'transcript_id': 'category',
                'genotype': 'category',
            })

        if library_id is not None:
            data['library_id'] = pd.Series([library_id], dtype="category")

        if sample_id is not None:
            data['sample_id'] = pd.Series([sample_id], dtype="category")

        snv_data.append(data)

    snv_data = scgenome.utils.concat_with_categories(snv_data, ignore_index=True)

    # Drop potential duplicates resulting from creating
    # sets of annotations from overlapping mutations
    # across different libraries
    snv_data = snv_data.drop_duplicates()

    return snv_data


def get_highest_snpeff_effect(snpeff_data):
    """ Select the highest ranked effect from snpeff data.
    """
    ordered_effect_impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

    ordered_effect_impacts = pd.DataFrame({
            'effect_impact': ordered_effect_impacts,
            'effect_impact_rank': range(len(ordered_effect_impacts))})
    ordered_effect_impacts['effect_impact'] = (
        ordered_effect_impacts['effect_impact'].astype(snpeff_data['effect_impact'].dtype))

    snpeff_data = snpeff_data.merge(ordered_effect_impacts)
    snpeff_data['coding_rank'] = (snpeff_data['amino_acid_change'].isnull() * 1)

    index_cols = ['chrom', 'coord', 'ref', 'alt']
    value_cols = ['gene_name', 'effect', 'effect_impact', 'amino_acid_change']

    snpeff_data = (
        snpeff_data[index_cols + value_cols + ['coding_rank', 'effect_impact_rank']]
        .sort_values(index_cols + ['coding_rank', 'effect_impact_rank'], ascending=True)
        .groupby(index_cols, sort=False, observed=True)
        .nth(0)
        .reset_index())

    snpeff_data = snpeff_data[[
        'chrom', 'coord', 'ref', 'alt',
        'gene_name', 'effect', 'effect_impact', 'amino_acid_change',
    ]]

    return snpeff_data


def load_snv_annotation_results(pseudobulk_dir, museq_filter=None, strelka_filter=None):
    """ Collate snv results into a single table.
    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    logging.info('starting load')
    mappability = load_snv_annotation_table(pseudobulk_dir, 'mappability')
    mappability = mappability[mappability['mappability'] > 0.99]

    strelka_results = load_snv_annotation_table(pseudobulk_dir, 'strelka')

    strelka_results = (
        strelka_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_strelka_score').reset_index())

    museq_results = load_snv_annotation_table(pseudobulk_dir, 'museq')
    museq_results = (
        museq_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_museq_score').reset_index())

    cosmic = load_snv_annotation_table(pseudobulk_dir, 'cosmic_status')
    logging.info(f'cosmic table with shape {cosmic.shape}, memory {cosmic.memory_usage().sum()}')
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = load_snv_annotation_table(pseudobulk_dir, 'snpeff')
    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info(f'snpeff table with shape {snpeff.shape}, memory {snpeff.memory_usage().sum()}')

    tnc = load_snv_annotation_table(pseudobulk_dir, 'trinuc')

    scgenome.utils.union_categories([
        mappability,
        cosmic,
        snpeff,
        tnc,
        strelka_results,
        museq_results,
    ])

    data = pd.concat([
        strelka_results.set_index(['chrom', 'coord', 'ref', 'alt']),
        museq_results.set_index(['chrom', 'coord', 'ref', 'alt']),
    ], axis=1).sort_index().reset_index()
    logging.info(f'merged snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(mappability, how='left')
    logging.info('post mappability with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(cosmic, how='left')
    logging.info('post cosmic with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(snpeff, how='left')
    logging.info('post snpeff with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(tnc, how='left')
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if museq_filter != -np.inf:
        data = data[data['max_museq_score'] > museq_filter]
        logging.info('post museq filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if strelka_filter != -np.inf:
        data = data[data['max_strelka_score'] > strelka_filter]
        logging.info('post strelka filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    logging.info('finishing load with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    for column in categorical_columns:
        data[column] = data[column].astype('category')

    logging.info(f'final snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    return data


def load_snv_data(
        results_dir,
        museq_filter=None,
        strelka_filter=None,
        positions=None,
        filter_sample_id=None,
        filter_library_id=None,
    ):
    """ Load filtered SNV annotation and count data
    
    Args:
        results_dir (str): results directory to load from.
        museq_score_threshold (float, optional): mutationseq score threshold. Defaults to None.
        strelka_score_threshold (float, optional): strelka score threshold. Defaults to None.

    Kwargs:
        filter_sample_id (str): restrict to specific sample id
        filter_library_id (str): restrict to specific library id

    Returns:
        pandas.DataFrame, pandas.DataFrame: SNV annotations, SNV counts
    """
    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['pseudobulk']

        snv_data = load_snv_annotation_results(
            pseudobulk_dir,
            museq_filter=museq_filter,
            strelka_filter=strelka_filter)

        assert not snv_data['coord'].isnull().any()

        positions = snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

        snv_count_data = load_snv_count_data(pseudobulk_dir, 'snv_union_counts.csv.gz', positions)
        snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']

        return {
            'snv_data': snv_data,
            'snv_count_data': snv_count_data,
        }

    if 'variant_calling' in analysis_dirs:
        variant_calling_dir = analysis_dirs['variant_calling']

        snv_data = load_snv_annotation_results(
            variant_calling_dir,
            museq_filter=museq_filter,
            strelka_filter=strelka_filter)

        assert not snv_data['coord'].isnull().any()

        return {
            'snv_data': snv_data,
        }

    if 'snv_genotyping' in analysis_dirs:
        variant_counting_dir = analysis_dirs['snv_genotyping']

        assert positions is not None

        snv_count_data = load_snv_count_data(
            variant_counting_dir, positions, 'counts.csv.gz',
            filter_sample_id=filter_sample_id,
            filter_library_id=filter_library_id)

        snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
        snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

        return {
            'snv_count_data': snv_count_data,
        }

    else:
        raise ValueError(f'no variant calling found for directory {results_dir}')

