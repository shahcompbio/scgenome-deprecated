import logging
import yaml
import os
import pandas as pd
import numpy as np

import scgenome.utils
import scgenome.loaders.utils


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


def load_snv_count_data(pseudobulk_dir, positions):
    """ Load per cell SNV count data
    
    Args:
        pseudobulk_dir (str): pseudobulk results directory
        positions (pandas.DataFrame): restrict to the specified positions
    
    Returns:
        pandas.DataFrame: SNV alt and ref counts per cell
    """
    snv_count_data = []

    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, 'snv_union_counts.csv.gz'):
        logging.info('Loading snv counts from {}'.format(filepath))

        data = pd.read_csv(
            filepath,
            dtype={
                'chrom': 'category',
                'ref': 'category',
                'alt': 'category',
                'cell_id': 'category',
            })

        logging.info(f'Loaded snv counts table with shape {data.shape}, size {data.memory_usage}')

        scgenome.utils.union_categories(
            [data, positions],
            cat_cols=['chrom', 'ref', 'alt'])
        data = data.merge(positions, how='inner')

        logging.info(f'Filtered snv counts table to shape {data.shape}, size {data.memory_usage}')

        snv_count_data.append(data)

    snv_count_data = scgenome.utils.concat_with_categories(snv_count_data, ignore_index=True)

    logging.info(f'Loaded all snv counts tables with shape {snv_count_data.shape}, size {snv_count_data.memory_usage}')

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
    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, f'snv_{table_name}.csv.gz'):
        logging.info(f'Loading snv {table_name} annotations from {filepath}')

        data = pd.read_csv(
            filepath,
            dtype={
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

        snv_data.append(data)

    snv_data = scgenome.utils.concat_with_categories(snv_data, ignore_index=True)

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

    index_cols = ['chrom', 'coord', 'ref', 'alt']
    value_cols = ['gene_name', 'effect', 'effect_impact', 'amino_acid_change']

    snpeff_data = (
        snpeff_data[index_cols + value_cols + ['effect_impact_rank']]
        .sort_values(index_cols + ['effect_impact_rank'], ascending=True)
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
    mappability['chrom'] = mappability['chrom'].astype(str)

    strelka_results = load_snv_annotation_table(pseudobulk_dir, 'strelka').rename(columns={'score': 'strelka_score'})
    for col in ('chrom', 'ref', 'alt'):
        strelka_results[col] = strelka_results[col].astype(str)

    museq_results = load_snv_annotation_table(pseudobulk_dir, 'museq').rename(columns={'score': 'museq_score'})
    for col in ('chrom', 'ref', 'alt'):
        museq_results[col] = museq_results[col].astype(str)

    cosmic = load_snv_annotation_table(pseudobulk_dir, 'cosmic_status')
    logging.info(f'cosmic table with shape {cosmic.shape}, size {cosmic.memory_usage}')
    cosmic['is_cosmic'] = True
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = load_snv_annotation_table(pseudobulk_dir, 'snpeff')
    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info(f'snpeff table with shape {snpeff.shape}, size {snpeff.memory_usage}')

    tnc = load_snv_annotation_table(pseudobulk_dir, 'trinuc')

    data = load_snv_annotation_table(pseudobulk_dir, 'allele_counts')
    logging.info(f'initial snv table with shape {data.shape}, size {data.memory_usage}')

    logging.info('summing snv counts')
    data = (
        data
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)[['alt_counts', 'ref_counts']]
        .sum().rename(columns={'alt_counts': 'alt_counts_sum', 'ref_counts': 'ref_counts_sum'}).reset_index())
    logging.info('total snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(mappability)
    logging.info('post mappability with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(cosmic, how='left')
    logging.info('post cosmic with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(snpeff, how='left')
    logging.info('post snpeff with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(tnc, how='left')
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(strelka_results, how='left')
    logging.info('post strelka with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    data = data.merge(museq_results, how='left')
    logging.info('post museq with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    if museq_filter != -np.inf:
        data = data[data['museq_score'] > museq_filter]
        logging.info('post museq filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    if strelka_filter != -np.inf:
        data = data[data['strelka_score'] > strelka_filter]
        logging.info('post strelka filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    logging.info('finishing load with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, size {data.memory_usage}')

    for column in categorical_columns:
        data[column] = data[column].astype('category')

    logging.info(f'final snv table with shape {data.shape}, size {data.memory_usage}')

    return data


def load_snv_data(
        results_dir,
        museq_filter=None,
        strelka_filter=None,
    ):
    """ Load filtered SNV annotation and count data
    
    Args:
        results_dir (str): results directory to load from.
        museq_score_threshold (float, optional): mutationseq score threshold. Defaults to None.
        strelka_score_threshold (float, optional): strelka score threshold. Defaults to None.
    
    Returns:
        pandas.DataFrame, pandas.DataFrame: SNV annotations, SNV counts
    """
    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'pseudobulk' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['pseudobulk']

    elif 'variant_calling' in analysis_dirs:
        pseudobulk_dir = analysis_dirs['variant_calling']

    else:
        raise ValueError(f'no pseudobulk found for directory {results_dir}')

    snv_data = load_snv_annotation_results(
        pseudobulk_dir,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter)

    assert not snv_data['coord'].isnull().any()

    positions = snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

    snv_count_data = load_snv_count_data(pseudobulk_dir, positions)
    snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
    snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

    return {
        'snv_data': snv_data,
        'snv_count_data': snv_count_data,
    }

