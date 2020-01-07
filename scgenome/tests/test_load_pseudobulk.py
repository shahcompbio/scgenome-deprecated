import os
import sys
import click
import logging
import collections

import dbclients.tantalus
import datamanagement.transfer_files

from scgenome.loaders.snv import load_snv_data
from scgenome.loaders.breakpoint import load_breakpoint_data
from scgenome.loaders.allele import load_haplotype_allele_data

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


dtypes_check = {
    'snv_data': {
        'chrom': 'category',
        'coord': 'int64',
        'ref': 'category',
        'alt': 'category',
        'mappability': 'float64',
        'is_cosmic': 'object',
        'gene_name': 'category',
        'effect': 'category',
        'effect_impact': 'category',
        'amino_acid_change': 'category',
        'tri_nucleotide_context': 'category',
        'max_strelka_score': 'float64',
        'max_museq_score': 'float64',
    },
    'snv_count_data': {
        'chrom': 'category',
        'coord': 'int64',
        'ref': 'category',
        'alt': 'category',
        'ref_counts': 'int64',
        'alt_counts': 'int64',
        'cell_id': 'category',
        'total_counts': 'int64',
    },
    'breakpoint_data': {
        'prediction_id': 'int64',
        'chromosome_1': 'object',
        'strand_1': 'object',
        'position_1': 'int64',
        'chromosome_2': 'object',
        'strand_2': 'object',
        'position_2': 'int64',
    },
    'breakpoint_count_data': {
        'prediction_id': 'int64',
        'cell_id': 'object',
        'read_count': 'int64',
    },
    'allele_counts': {
        'allele_id': 'int64',
        'cell_id': 'category',
        'chromosome': 'category',
        'end': 'int64',
        'hap_label': 'int64',
        'readcount': 'int64',
        'start': 'int64',
    },
}


def test_pseudobulk_data(snv_results_tables, breakpoint_results_tables, haplotype_results_tables):
    failed = False
    for results_tables in (snv_results_tables, breakpoint_results_tables, haplotype_results_tables):
        for table_name, table_data in results_tables.items():
            logging.info(f'table {table_name} has size {len(table_data)}')
            for column_name, dtype_name in dtypes_check[table_name].items():
                column_dtype = str(results_tables[table_name][column_name].dtype)
                expected_dtype_names = (dtype_name,)
                if dtype_name == 'int64':
                    expected_dtype_names = ('int64', 'Int64')
                if not column_dtype in expected_dtype_names:
                    logging.error(f'{column_name} has dtype {column_dtype} not any of {expected_dtype_names}')
                    failed = True
    if failed:
        raise Exception('failed dtype check')


def test_load_local_pseudobulk_data(results_dir):
    snv_results_tables = load_snv_data(
        results_dir,
    )

    breakpoint_results_tables = load_breakpoint_data(
        results_dir,
    )

    haplotype_results_tables = load_haplotype_allele_data(
        results_dir,
    )

    test_pseudobulk_data(
        snv_results_tables,
        breakpoint_results_tables,
        haplotype_results_tables,
    )

    logging.info(f'successfully loaded results from {results_dir}')


def test_load_stored_pseudobulk_data(tantalus_api, ticket_id, local_cache_directory=None, local_storage_name=None):
    if local_cache_directory is not None and local_storage_name is not None:
        raise ValueError('local_cache_directory and local_storage_name are mutually exclusive')

    if local_cache_directory is None and local_storage_name is None:
        raise ValueError('require one of local_cache_directory and local_storage_name')

    if local_cache_directory is not None:
        ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

        for results in ticket_results:
            filepaths = datamanagement.transfer_files.cache_dataset(
                tantalus_api,
                results['id'],
                'resultsdataset',
                'singlecellresults',
                local_cache_directory,
            )

        local_results_directory = local_cache_directory

    elif local_storage_name is not None:
        local_results_directory = tantalus_api.get('storage', name=local_storage_name)['storage_directory']

    ticket_directory = os.path.join(local_results_directory, ticket_id)

    test_load_local_pseudobulk_data(ticket_directory)


@click.group()
def cli():
    pass


@cli.command()
@click.argument('ticket_id')
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
def test_cached_single_ticket(ticket_id, local_cache_directory=None, local_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    test_load_stored_pseudobulk_data(
        tantalus_api,
        ticket_id,
        local_cache_directory=local_cache_directory,
        local_storage_name=local_storage_name,
    )


@cli.command()
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
def test_cached_multi_ticket(ticket_id, local_cache_directory=None, local_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    hmmcopy_analyses = tantalus_api.list('analysis', analysis_type__name='hmmcopy')

    version_tickets = collections.defaultdict(list)
    for analysis in hmmcopy_analyses:
        version_tickets[analysis['version']].append(analysis['jira_ticket'])
    
    ticket_ids = []
    for version in version_tickets:
        ticket_ids.append(version_tickets[version][-1])

    for ticket_id in ticket_ids:
        logging.info(ticket_id)
        test_load_stored_pseudobulk_data(
            tantalus_api,
            ticket_id,
            local_cache_directory=local_cache_directory,
            local_storage_name=local_storage_name,
        )


@cli.command()
@click.argument('results_directory')
def test_local_results(results_directory):
    test_load_local_pseudobulk_data(results_directory)


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    cli()

