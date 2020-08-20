import logging
import yaml
import os
import packaging
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

def load_snv_count_data_from_filenames(files, positions, filter_sample_id=None, 
    filter_library_id=None):
    return _process_snv_count_data(scgenome.loaders.utils._prep_filenames_for_loading(files),
        positions, filter_sample_id=filter_sample_id, filter_library_id=filter_library_id
    )


def load_snv_count_data(pseudobulk_dir, suffix, positions, filter_sample_id=None, 
    filter_library_id=None):
    """ Load per cell SNV count data
    
    Args:
        pseudobulk_dir (str): pseudobulk results directory
        suffix (str): suffix of snv count tables
        positions (pandas.DataFrame): restrict to the specified positions

    Kwargs:
        filter_sample_id (str): restrict to specific sample id
        filter_library_id (str): restrict to specific library id
        files: (list of str): optionally pass list of counts filepaths too selectively load count data
    Returns:
        pandas.DataFrame: SNV alt and ref counts per cell
    """

    files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, suffix)

    return _process_snv_count_data(files, positions, filter_sample_id=filter_sample_id, 
        filter_library_id=filter_library_id
    )


def _process_snv_count_data(files, positions, filter_sample_id=None, filter_library_id=None):

    snv_count_data = []

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

    logging.info(f'Loaded all snv counts tables with shape {snv_count_data.shape}, memory \
        {snv_count_data.memory_usage().sum()}')

    return snv_count_data


def load_snv_annotation_table(files):
    """ Load SNV annotation data

    Args:
        pseudobulk_dir (str): pseudobulk results directory
        table_name (str): name of annotation table to load.

    Returns:
        pandas.DataFrame: SNVs annotation data per sample / library
    """
    snv_data = []

    for sample_id, library_id, filepath in files:
        logging.info(f'Loading from {filepath}')

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


def load_snv_annotation_results_from_filenames(mappability_path, strelka_path, museq_path, cosmic_path, 
    snpeff_path, dbsnp_path, trinuc_path, museq_filter=None, strelka_filter=None
):
    """ Collate snv results into a single table from input filenames. path inputs must be lists of file strings. 
    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    logging.info('starting load')


    mappability = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(
            mappability_path))
    mappability = mappability[mappability['mappability'] > 0.99]

    strelka_results = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(
            strelka_path))

    strelka_results = (
        strelka_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_strelka_score').reset_index())

    museq_results = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(
            museq_path))
    museq_results = (
        museq_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_museq_score').reset_index())

    cosmic = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(cosmic_path))
    logging.info(f'cosmic table with shape {cosmic.shape}, memory {cosmic.memory_usage().sum()}')
    cosmic['is_cosmic'] = 1
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(snpeff_path))

    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info(f'snpeff table with shape {snpeff.shape}, memory {snpeff.memory_usage().sum()}')

    dbsnp = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(dbsnp_path))
    dbsnp = dbsnp[['chrom', 'coord', 'ref', 'alt', 'exact_match']].rename(columns={'exact_match': 'is_dbsnp'})
    logging.info(f'dbsnp table with shape {dbsnp.shape}, memory {dbsnp.memory_usage().sum()}')

    tnc = load_snv_annotation_table(scgenome.loaders.utils._prep_filenames_for_loading(trinuc_path))

    return _concat_annotation_results(mappability, cosmic, snpeff, dbsnp, tnc, strelka_results, 
        museq_results, museq_filter=museq_filter, strelka_filter=strelka_filter
    )


def _concat_annotation_results(mappability, cosmic, snpeff, dbsnp, tnc, strelka_results, museq_results, 
    museq_filter=None, strelka_filter=None
):
    '''
    private function to concatenate and filter snv annotation data
    '''

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


def load_snv_annotation_results(pseudobulk_dir, museq_filter=None, strelka_filter=None):
    """ Collate snv results into a single table.
    """

    if museq_filter is None:
        museq_filter = default_museq_filter

    if strelka_filter is None:
        strelka_filter = default_strelka_filter

    logging.info('starting load')

    mappability_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_mappability.csv.gz'
    )
    mappability = load_snv_annotation_table(mappability_files)
    mappability = mappability[mappability['mappability'] > 0.99]

    strelka_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_strelka.csv.gz'
    )
    strelka_results = load_snv_annotation_table(strelka_files)

    strelka_results = (
        strelka_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_strelka_score').reset_index())

    museq_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_museq.csv.gz'
    )
    museq_results = load_snv_annotation_table(museq_files)
    museq_results = (
        museq_results
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['score']
        .max().rename('max_museq_score').reset_index())

    cosmic_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_cosmic_status.csv.gz'
    )
    cosmic = load_snv_annotation_table(cosmic_files)
    logging.info(f'cosmic table with shape {cosmic.shape}, memory {cosmic.memory_usage().sum()}')
    cosmic['is_cosmic'] = 1
    cosmic = cosmic[['chrom', 'coord', 'ref', 'alt', 'is_cosmic']].drop_duplicates()

    snpeff_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_snpeff.csv.gz'
    )

    snpeff = load_snv_annotation_table(snpeff_files)

    snpeff = get_highest_snpeff_effect(snpeff)
    logging.info(f'snpeff table with shape {snpeff.shape}, memory {snpeff.memory_usage().sum()}')

    dbsnp_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_dbsnp_status.csv.gz'
    )
    dbsnp = load_snv_annotation_table(dbsnp_files)
    dbsnp = dbsnp[['chrom', 'coord', 'ref', 'alt', 'exact_match']].rename(columns={'exact_match': 'is_dbsnp'})
    logging.info(f'dbsnp table with shape {dbsnp.shape}, memory {dbsnp.memory_usage().sum()}')

    trinuc_files = scgenome.loaders.utils.get_pseudobulk_files(
        pseudobulk_dir, 'snv_trinuc.csv.gz'
    )
    tnc = load_snv_annotation_table(trinuc_files)

    return _concat_annotation_results(mappability, cosmic, snpeff, dbsnp, tnc, strelka_results, museq_results, 
        museq_filter=museq_filter, strelka_filter=strelka_filter
    )


def load_snv_data_from_files(        
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
    positions=None,
    filter_sample_id=None,
    filter_library_id=None,
    snv_annotation=False,
    snv_counts=False
):

    """ Load filtered SNV annotation and count data
    """

    outputs = {}

    if snv_annotation:
        snv_data = load_snv_annotation_results_from_filenames(
            mappability_path, strelka_path, museq_path, cosmic_path, snpeff_path, dbsnp_path, trinuc_path, 
            museq_filter=museq_filter, strelka_filter=strelka_filter
        )

        assert not snv_data['coord'].isnull().any()

        outputs["snv_data"] = snv_data

    if snv_counts:
        assert positions is not None

        snv_count_data = load_snv_count_data_from_filenames(
            counts_path,
            positions,
            filter_sample_id=filter_sample_id,
            filter_library_id=filter_library_id)

        snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
        snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

        outputs["snv_count_data"] = snv_count_data

    return outputs


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

        if len(variant_calling_dir) == 0:
            raise ValueError(f'found {len(variant_calling_dir)} dirs for variant_calling')

        elif len(variant_calling_dir) > 1:
            if filter_sample_id is None:
                raise ValueError(f'found {len(variant_calling_dir)} without filter_sample_id')

            filtered_variant_calling_dir = list(filter(lambda a: f'sample_{filter_sample_id}' in a, variant_calling_dir))

            if len(filtered_variant_calling_dir) != 1:
                raise ValueError(f'found {len(filtered_variant_calling_dir)} in {variant_calling_dir} matching filter_sample_id')

            variant_calling_dir = filtered_variant_calling_dir

        snv_data = load_snv_annotation_results(
            variant_calling_dir[0],
            museq_filter=museq_filter,
            strelka_filter=strelka_filter)

        assert not snv_data['coord'].isnull().any()

        return {
            'snv_data': snv_data,
        }

    if 'snv_genotyping' in analysis_dirs:
        assert positions is not None

        variant_counting_dir = analysis_dirs['snv_genotyping']

        if len(variant_counting_dir) != 1:
            raise ValueError(f'found {len(variant_counting_dir)} dirs for snv_genotyping')
        variant_counting_dir = variant_counting_dir[0]

        manifest_filename = os.path.join(variant_counting_dir, 'metadata.yaml')
        manifest = yaml.load(open(manifest_filename))

        suffix = 'counts.csv.gz'
        if packaging.version.parse(manifest['meta']['version']) > packaging.version.parse('v0.6.0'):
            if filter_sample_id is None or filter_library_id is None:
                raise ValueError('both filter_sample_id and filter_library_id must be specified')
            suffix = f'{filter_sample_id}_{filter_library_id}_counts.csv.gz'

        snv_count_data = load_snv_count_data(
            variant_counting_dir,
            suffix,
            positions,
            filter_sample_id=filter_sample_id,
            filter_library_id=filter_library_id)

        snv_count_data['total_counts'] = snv_count_data['ref_counts'] + snv_count_data['alt_counts']
        snv_count_data['sample_id'] = snv_count_data['cell_id'].apply(lambda a: a.split('-')[0]).astype('category')

        return {
            'snv_count_data': snv_count_data,
        }

    else:
        raise ValueError(f'no variant calling found for directory {results_dir}')

