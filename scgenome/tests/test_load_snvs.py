import sys
import click
import logging
import collections

import dbclients.tantalus
import datamanagement.transfer_files

from scgenome.loaders.snv import load_cached_snv_data

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


@click.command()
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
@click.option('--single_ticket_id')
def test_load_cached_snv_data(local_cache_directory, local_storage_name, single_ticket_id):
    if local_cache_directory is not None and local_storage_name is not None:
        raise ValueError('local_cache_directory and local_storage_name are mutually exclusive')

    if local_cache_directory is None and local_storage_name is None:
        raise ValueError('require one of local_cache_directory and local_storage_name')

    tantalus_api = dbclients.tantalus.TantalusApi()

    if single_ticket_id is None:
        pseudobulk_analyses = tantalus_api.list('analysis', analysis_type__name='pseudobulk')

        version_tickets = collections.defaultdict(list)
        for analysis in pseudobulk_analyses:
            version_tickets[analysis['version']].append(analysis['jira_ticket'])
        
        ticket_ids = []
        for version in version_tickets:
            ticket_ids.append(version_tickets[version][-1])

    else:
        ticket_ids = (single_ticket_id,)

    for ticket_id in ticket_ids:
        logging.info(ticket_id)

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

        results_tables = load_cached_snv_data(
            ticket_id,
            local_results_directory,
        )

        for table_name, table_data in results_tables.items():
            print(table_name)
            print(table_data.shape)
            print(table_data.head())


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    test_load_cached_snv_data()
