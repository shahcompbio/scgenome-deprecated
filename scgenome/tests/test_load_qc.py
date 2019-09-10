import sys
import click
import logging
import collections

import dbclients.tantalus

from scgenome.loaders.qc import load_cached_qc_data
from scgenome.loaders.align import load_cached_align_data
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
        'end': 'float64',
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
    },
    'gc_metrics': {
    },
}


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
            results_storage_name='singlecellresults',
        )

        results_tables = load_cached_qc_data(
            ticket_id,
            local_cache_directory,
        )

        for table_name, table_data in results_tables.items():
            logging.info(f'table {table_name} has size {len(table_data)}')
            for column_name, dtype_name in dtypes_check[table_name].items():
                column_dtype = str(results_tables[table_name][column_name].dtype)
                if not column_dtype == dtype_name:
                    raise Exception(f'{column_name} has dtype {column_dtype} not {dtype_name}')


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    test_load_cached_qc_data()
