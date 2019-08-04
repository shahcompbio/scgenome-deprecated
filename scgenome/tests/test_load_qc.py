import sys
import click
import logging
import collections

import dbclients.tantalus

from scgenome.loaders.qc import load_cached_qc_data
from scgenome.loaders.align import load_cached_align_data
from scgenome.db.qc import cache_qc_results


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


@click.command()
@click.argument('local_cache_directory')
@click.option('--single_ticket_id')
def test_load_cached_qc_data(local_cache_directory, single_ticket_id):
    if single_ticket_id is None:
        tantalus_api = dbclients.tantalus.TantalusApi()

        hmmcopy_analyses = tantalus_api.list('analysis', analysis_type__name='hmmcopy')

        version_tickets = collections.defaultdict(list)
        for analysis in hmmcopy_analyses:
            version_tickets[analysis['version']].append(analysis['jira_ticket'])
        
        ticket_ids = []
        for version in version_tickets:
            ticket_ids.append(version_tickets[version][-1])

    else:
        ticket_ids = (single_ticket_id,)

    for ticket_id in ticket_ids:
        logging.info(ticket_id)

        cache_qc_results(
            ticket_id,
            local_cache_directory,
        )

        results_tables = load_cached_qc_data(
            ticket_id,
            local_cache_directory,
        )

        for table_name, table_data in results_tables.items():
            print(table_name)
            print(table_data.shape)
            print(table_data.head())


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    test_load_cached_qc_data()
