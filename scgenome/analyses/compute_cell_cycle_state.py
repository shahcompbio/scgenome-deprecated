import os
import sys
import click
import logging
import shutil
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from datamanagement.utils.utils import make_dirs
import scgenome.dataimport
import scgenome.cellcycle


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


analysis_type = 'cell_state_classifier'
analysis_version = 'v0.0.1'


results_type = 'cell_state_prediction'
results_version = 'v0.0.1'


def get_unprocessed_hmmcopy(tantalus_api):
    unprocessed = {}

    for results in tantalus_api.list('results', results_type='hmmcopy'):
        # Retrieve the jira ticket of the hmmcopy analysis
        hmmcopy_analysis = tantalus_api.get('analysis', id=results['analysis'])
        jira_ticket = hmmcopy_analysis['jira_ticket']

        # Check for an existing analysis with this hmmcopy as input
        try:
            analysis = tantalus_api.get(
                'analysis',
                analysis_type=analysis_type,
                version=analysis_version,
                input_results__id=results['id'],
            )
        except dbclients.basicclient.NotFoundError:
            analysis = None

        if analysis is not None:
            logging.info('hmmcopy ticket {} has cell cycle analysis {}'.format(
                jira_ticket, analysis['name']))
            continue

        assert jira_ticket not in unprocessed
        unprocessed[jira_ticket] = results

    return unprocessed


def run_analysis(
        tantalus_api, hmmcopy_results, jira_ticket, hmmcopy_jira_ticket,
        cache_directory, results_storage_name, archive_storage_name=None):

    assert len(hmmcopy_results['libraries']) == 1
    library_id = hmmcopy_results['libraries'][0]['library_id']
    library_pk = hmmcopy_results['libraries'][0]['id']

    results_storage = tantalus_api.get(
        'storage',
        name=results_storage_name,
    )

    results_filename = os.path.join(
        'single_cell_indexing',
        '{results_type}',
        '{results_version}',
        '{jira_ticket}',
        '{hmmcopy_jira_ticket}_{library_id}.csv',
    ).format(
        results_type=results_type,
        results_version=results_version,
        hmmcopy_jira_ticket=hmmcopy_jira_ticket,
        jira_ticket=jira_ticket,
        library_id=library_id,
    )

    results_filepath = os.path.join(
        results_storage['storage_directory'],
        results_filename,
    )

    make_dirs(os.path.dirname(results_filepath))

    logging.info('loading data for hmmcopy ticket {}'.format(
        hmmcopy_jira_ticket))

    results = scgenome.dataimport.import_cn_data(
        [hmmcopy_jira_ticket],
        cache_directory,
    )

    cn_data = results['hmmcopy_reads']

    logging.info('calculating cell cycle state')

    cell_cycle_data = scgenome.cellcycle.predict(cn_data)
    cell_cycle_data.to_csv(results_filepath, index=False)

    logging.info('registering results with tantalus')

    analysis_name = '{}_{}_{}_{}'.format(
        analysis_type, analysis_version,
        hmmcopy_jira_ticket, library_id,
    )

    analysis = tantalus_api.get_or_create(
        'analysis',
        name=analysis_name,
        analysis_type=analysis_type,
        version=analysis_version,
        jira_ticket=jira_ticket,
        status='complete',
        args={},
        input_results=[hmmcopy_results['id']],
    )

    results_name = '{}_{}_{}_{}'.format(
        results_type, results_version,
        hmmcopy_jira_ticket, library_id,
    )

    results_file_resource, results_file_instance = tantalus_api.add_file(
        results_storage_name, results_filepath, update=True)
    results_file_pk = results_file_resource['id']

    results = tantalus_api.get_or_create(
        'results',
        name=results_name,
        results_type=results_type,
        results_version=results_version,
        libraries=[library_pk],
        analysis=analysis['id'],
        file_resources=[results_file_pk],
    )

    if archive_storage_name is not None:
        datamanagement.transfer_files.transfer_dataset(
            tantalus_api,
            results['id'],
            'resultsdataset',
            results_storage_name,
            archive_storage_name,
        )


@click.command()
@click.argument('jira_ticket', nargs=1)
@click.argument('cache_directory', nargs=1)
@click.argument('results_storage_name', nargs=1)
@click.option('--archive_storage_name', required=False)
def run_all_analyses(jira_ticket, cache_directory, results_storage_name, archive_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    datasets = get_unprocessed_hmmcopy(tantalus_api)
    for hmmcopy_jira_ticket, dataset in datasets.items():
        logging.info('processing dataset {}'.format(dataset['name']))

        try:
            run_analysis(
                tantalus_api, dataset, jira_ticket, hmmcopy_jira_ticket,
                cache_directory, results_storage_name,
                archive_storage_name=archive_storage_name)
        except Exception as e:
            logging.warning('processing of dataset {} failed with exception {}'.format(dataset['name'], e))

        ticket_directory = os.path.join(cache_directory, hmmcopy_jira_ticket)
        logging.info('cleaning ticket directory {}'.format(ticket_directory))
        try:
            shutil.rmtree(ticket_directory)
        except Exception as e:
            logging.warning('unable to remove {}, exception {}'.format(ticket_directory, e))


if __name__ == "__main__":
    # Set up the root logger
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    run_all_analyses()

