import logging
import os

import numpy as np
import packaging
import pandas as pd
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils
import yaml

default_museq_filter = 0.9
default_strelka_filter = 20.
default_mappability_filter = 0.99

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


def load_snv_genotyping_results(results_dir, positions=None):
    """ Load per cell SNV count data from results directory
    
    Args:
        results_dir (str): results directory
        positions (pandas.DataFrame): restrict to the specified positions

    Returns:
        pandas.DataFrame: SNV alt and ref counts per cell
    """

    counts_filepath = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'counts.csv.gz', analysis_type='snv_genotyping')

    return load_snv_genotyping_files(counts_filepath, positions=positions)


def load_snv_genotyping_files(filepath, positions=None):
    """ Load per cell SNV count data from filepath
    
    Args:
        filepath (str): counts filepath
        positions (pandas.DataFrame): restrict to the specified positions

    Returns:
        pandas.DataFrame: SNV alt and ref counts per cell
    """

    logging.info('Loading snv counts from {}'.format(filepath))

    data = []

    csv_input = scgenome.csvutils.CsvInput(filepath)

    chunk_iter = csv_input.read_csv(
        chunksize=10 ** 6,
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
        if positions is not None:
            scgenome.utils.union_categories(
                [chunk, positions],
                cat_cols=['chrom', 'ref', 'alt'])
            chunk = chunk.merge(positions)

        data.append(chunk)

    data = scgenome.utils.concat_with_categories(data, ignore_index=True)

    logging.info(f'Loaded snv counts table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if positions is not None:
        scgenome.utils.union_categories(
            [data, positions],
            cat_cols=['chrom', 'ref', 'alt'])
        data = data.merge(positions, how='inner')

    logging.info(f'Filtered snv counts table to shape {data.shape}, memory {data.memory_usage().sum()}')

    shape = data.shape
    memory = data.memory_usage().sum()
    logging.info(f'Loaded all snv counts tables with shape {shape}, memory {memory}')

    data['total_counts'] = data['ref_counts'] + data['alt_counts']

    return {'snv_count_data': data}


def load_snv_annotation_table(filepath):
    """ Load SNV annotation data

    Args:
        filepath (str): path of snv annotation table.

    Returns:
        pandas.DataFrame: SNVs annotation data
    """
    logging.info(f'Loading from {filepath}')

    csv_input = scgenome.csvutils.CsvInput(filepath)
    snv_data = csv_input.read_csv(
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

    snpeff_data = (snpeff_data[index_cols + value_cols + ['coding_rank', 'effect_impact_rank']]
                   .sort_values(index_cols + ['coding_rank', 'effect_impact_rank'], ascending=True)
                   .groupby(index_cols, sort=False, observed=True)
                   .nth(0)
                   .reset_index()
                   )

    snpeff_data = snpeff_data[[
        'chrom', 'coord', 'ref', 'alt',
        'gene_name', 'effect', 'effect_impact', 'amino_acid_change',
    ]]

    return snpeff_data


def load_snv_annotation_results(
        results_dir,
        museq_filter=None,
        strelka_filter=None,
        mappability_filter=None,
    ):
    """ Load per cell SNV count data from results directory
    
    Args:
        results_dir (str): results directory

    Kwargs:
        museq_filter (float, optional): mutationseq score threshold. Defaults to None.
        strelka_filter (float, optional): strelka score threshold. Defaults to None.
        mappability_filter (float, optional): mappability threshold. Defaults to None.

    Returns:
        pandas.DataFrame: SNV alt and ref per cell
    """

    mappability_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_mappability.csv.gz', analysis_type='variant_calling')

    strelka_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_strelka.csv.gz', analysis_type='variant_calling')

    museq_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_museq.csv.gz', analysis_type='variant_calling')

    cosmic_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_cosmic_status.csv.gz', analysis_type='variant_calling')

    snpeff_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_snpeff.csv.gz', analysis_type='variant_calling')

    dbsnp_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_dbsnp_status.csv.gz', analysis_type='variant_calling')

    trinuc_path = scgenome.loaders.utils.find_results_filepath(
        results_dir, 'snv_trinuc.csv.gz', analysis_type='variant_calling')

    return load_snv_annotation_files(
        mappability_path,
        strelka_path,
        museq_path,
        cosmic_path,
        snpeff_path,
        dbsnp_path,
        trinuc_path,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter,
        mappability_filter=mappability_filter,
    )


def load_snv_annotation_files(
        mappability_path,
        strelka_path,
        museq_path,
        cosmic_path,
        snpeff_path,
        dbsnp_path,
        trinuc_path,
        museq_filter=None,
        strelka_filter=None,
        mappability_filter=None,
    ):
    """ Collate snv results into a single table from input filenames. path inputs must be lists of file strings. 

    Args:
        mappability_path (str): mappability file path
        strelka_path (str): strelka file path
        museq_path (str): museq file path
        cosmic_path (str): cosmic file path
        snpeff_path (str): snpeff file path
        dbsnp_path (str): dbsnp file path
        trinuc_path (str): trinuc file path

    Kwargs:
        museq_filter (float, optional): mutationseq score threshold. Defaults to None.
        strelka_filter (float, optional): strelka score threshold. Defaults to None.
        mappability_filter (float, optional): mappability threshold. Defaults to None.

    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    if mappability_filter is None:
        mappability_filter = default_mappability_filter

    logging.info('starting load')

    mappability = load_snv_annotation_table(mappability_path)

    strelka_results = load_snv_annotation_table(strelka_path)
    strelka_results = (
        strelka_results
            .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
            .max().rename('max_strelka_score').reset_index())

    museq_results = load_snv_annotation_table(museq_path)
    museq_results = (
        museq_results
            .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
            .max().rename('max_museq_score').reset_index())

    cosmic = load_snv_annotation_table(cosmic_path)
    logging.info(f'cosmic table with shape {cosmic.shape}, memory {cosmic.memory_usage().sum()}')
    cosmic['is_cosmic'] = 1
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = load_snv_annotation_table(snpeff_path)

    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info(f'snpeff table with shape {snpeff.shape}, memory {snpeff.memory_usage().sum()}')

    dbsnp = load_snv_annotation_table(dbsnp_path)
    dbsnp = dbsnp.query('exact_match == 1')[['chrom', 'coord', 'ref', 'alt']].drop_duplicates().assign(is_dbsnp=1)
    logging.info(f'dbsnp table with shape {dbsnp.shape}, memory {dbsnp.memory_usage().sum()}')

    tnc = load_snv_annotation_table(trinuc_path)

    scgenome.utils.union_categories([
        mappability,
        cosmic,
        snpeff,
        dbsnp,
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
    data['is_cosmic'] = data['is_cosmic'].fillna(0).astype(int)
    logging.info('post cosmic with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(snpeff, how='left')
    logging.info('post snpeff with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(dbsnp, how='left')
    data['is_dbsnp'] = data['is_dbsnp'].fillna(0).astype(int)
    logging.info('post dbsnp with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    data = data.merge(tnc, how='left')
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if museq_filter is not None:
        data = data[data['max_museq_score'] >= museq_filter]
        logging.info('post museq filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if strelka_filter is not None:
        data = data[data['max_strelka_score'] >= strelka_filter]
        logging.info(
            'post strelka filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    if mappability_filter is not None:
        data = data[data['mappability'] >= mappability_filter]
        logging.info(
            'post mappability filter with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
        logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    logging.info('finishing load with snv count {}'.format(data[['chrom', 'coord']].drop_duplicates().shape[0]))
    logging.info(f'snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    for column in categorical_columns:
        data[column] = data[column].astype('category')

    logging.info(f'final snv table with shape {data.shape}, memory {data.memory_usage().sum()}')

    return {'snv_data': data}


def load_snv_results(
        variant_calling_results_dir,
        snv_genotyping_results_dir,
        museq_filter=None,
        strelka_filter=None,
        mappability_filter=None,
    ):
    """ Load filtered SNV annotation and genotyping data
    """

    variant_calling_results_dir = scgenome.loaders.utils.find_results_directory(
        variant_calling_results_dir, 'variant_calling')

    results_tables = load_snv_annotation_results(
        variant_calling_results_dir,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter,
        mappability_filter=mappability_filter,
    )

    snv_genotyping_results_dir = scgenome.loaders.utils.find_results_directory(
        snv_genotyping_results_dir, 'snv_genotyping')

    positions = results_tables['snv_data'][['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

    snv_genotyping_results = load_snv_genotyping_results(snv_genotyping_results_dir, positions)
    results_tables.update(snv_genotyping_results)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_snv_files(
        mappability_path,
        strelka_path,
        museq_path,
        cosmic_path,
        snpeff_path,
        dbsnp_path,
        trinuc_path,
        counts_path,
        museq_filter=None,
        strelka_filter=None,
        mappability_filter=None,
    ):
    """ Load filtered SNV annotation and genotyping data
    """

    results_tables = load_snv_annotation_files(
        mappability_path,
        strelka_path,
        museq_path,
        cosmic_path,
        snpeff_path,
        dbsnp_path,
        trinuc_path,
        museq_filter=museq_filter,
        strelka_filter=strelka_filter,
        mappability_filter=mappability_filter,
    )

    assert not results_tables['snv_data']['coord'].isnull().any()

    positions = results_tables['snv_data'][['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

    snv_genotyping_results = load_snv_genotyping_files(counts_path, positions)
    results_tables.update(snv_genotyping_results)

    return results_tables

