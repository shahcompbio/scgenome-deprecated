import os
import sys
import click
import logging
import collections

import dbclients.tantalus
import datamanagement.transfer_files

from scgenome.loaders.qc import load_qc_data
from scgenome.db.qc import cache_qc_results


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


dtypes_check = {
    'align_metrics': {
        'cell_id': 'category',
        'total_mapped_reads': 'int64',
        'library_id': 'category',
    },
    'hmmcopy_reads': {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'reads': 'int64',
        'gc': 'float64',
        'copy': 'float64',
        'state': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_segs': {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'state': 'int64',
        'median': 'float64',
        'multiplier': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_metrics': {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'order': 'float64',
        'total_mapped_reads_hmmcopy': 'int64',
        'experimental_condition': 'object',
        'mean_copy': 'float64',
        'state_mode': 'int64',
    },
    'annotation_metrics': {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'total_mapped_reads_hmmcopy': 'int64',
        'mean_copy': 'float64',
        'state_mode': 'int64',
        'order': 'float64',
        'experimental_condition': 'object',
        'quality': 'float64',
        'is_s_phase': 'boolean',
    },
    'gc_metrics': {
    },
}


def test_qc_data(results_tables):
    for table_name, table_data in results_tables.items():
        logging.info(f'table {table_name} has size {len(table_data)}')
        for column_name, dtype_name in dtypes_check[table_name].items():
            column_dtype = str(results_tables[table_name][column_name].dtype)
            if not column_dtype == dtype_name:
                raise Exception(f'{column_name} has dtype {column_dtype} not {dtype_name}')


def test_load_local_qc_data(results_dir):
    results_tables = load_qc_data(results_dir)
    test_qc_data(results_tables)


def test_load_stored_qc_data(tantalus_api, ticket_id, local_cache_directory=None, local_storage_name=None):
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

    test_load_local_qc_data(ticket_directory)


@click.group()
def cli():
    pass


@cli.command()
@click.argument('ticket_id')
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
def test_cached_single_ticket(ticket_id, local_cache_directory=None, local_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    test_load_stored_qc_data(
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
        test_load_stored_qc_data(
            tantalus_api,
            ticket_id,
            local_cache_directory=local_cache_directory,
            local_storage_name=local_storage_name,
        )


@cli.command()
@click.argument('results_directory')
def test_local_results(results_directory):
    test_load_local_qc_data(results_directory)


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    cli()
