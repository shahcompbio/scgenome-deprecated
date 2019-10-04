import os
import sys
import click
import logging
import shutil
import yaml
import pandas as pd

import cell_cycle_classifier.api
import dbclients.tantalus
import datamanagement.transfer_files
from datamanagement.utils.utils import make_dirs

from scgenome.loaders.align import load_align_data
from scgenome.loaders.hmmcopy import load_hmmcopy_data


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


analysis_type = 'cell_state_classifier'
analysis_version = 'v0.0.3'


results_type = 'cell_state_prediction'
results_version = 'v0.0.3'


def get_unprocessed_hmmcopy(tantalus_api, hmmcopy_tickets=None, rerun=False):
    if hmmcopy_tickets is None or len(hmmcopy_tickets) == 0:
        hmmcopy_results = list(tantalus_api.list('results', results_type='hmmcopy'))

    else:
        hmmcopy_results = []
        for ticket in hmmcopy_tickets:
            hmmcopy_results.append(tantalus_api.get('results', results_type='hmmcopy', analysis__jira_ticket=ticket))

    unprocessed = {}

    for results in hmmcopy_results:
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

        if analysis is not None and not rerun:
            logging.info('hmmcopy ticket {} has cell cycle analysis {}'.format(
                jira_ticket, analysis['name']))
            continue

        assert jira_ticket not in unprocessed
        unprocessed[jira_ticket] = results

    return unprocessed


def run_analysis(
        tantalus_api, hmmcopy_results, jira_ticket,
        results_storage_name, archive_storage_name=None, rerun=False):

    assert len(hmmcopy_results['libraries']) == 1
    library_id = hmmcopy_results['libraries'][0]['library_id']
    library_pk = hmmcopy_results['libraries'][0]['id']

    results_storage = tantalus_api.get(
        'storage',
        name=results_storage_name,
    )

    results_dir = os.path.join(
        '{jira_ticket}',
        'results',
        '{results_type}',
    ).format(
        jira_ticket=jira_ticket,
        results_type=results_type,
    )

    relative_filename = f'{library_id}_{results_type}.csv'

    results_filename = os.path.join(
        results_dir,
        relative_filename,
    )

    results_filepath = os.path.join(
        results_storage['storage_directory'],
        results_filename,
    )

    metadata_filename = os.path.join(
        results_dir,
        'metadata.yaml',
    )

    metadata_filepath = os.path.join(
        results_storage['storage_directory'],
        metadata_filename,
    )

    make_dirs(os.path.dirname(results_filepath))

    logging.info('loading data for hmmcopy ticket {}'.format(
        jira_ticket))

    hmmcopy_results_dir = os.path.join(results_storage['storage_directory'], jira_ticket)

    results = load_align_data(hmmcopy_results_dir)
    results.update(load_hmmcopy_data(hmmcopy_results_dir))

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']
    align_metrics_data = results['align_metrics']

    logging.info('calculating cell cycle state')

    cell_cycle_data = cell_cycle_classifier.api.train_classify(cn_data, metrics_data, align_metrics_data)
    cell_cycle_data.to_csv(results_filepath, index=False)

    metadata = {
        'filenames': [
            relative_filename,
        ],
        'meta': {
            'library_id': hmmcopy_results['libraries'][0]['library_id'],
            'sample_id': [a['sample_id'] for a in hmmcopy_results['samples']],
            'type': results_type,
            'version': results_version,
        }
    }

    with open(metadata_filepath, 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)

    logging.info('registering results with tantalus')

    analysis_name = '{}_{}_{}_{}'.format(
        analysis_type, analysis_version,
        jira_ticket, library_id,
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
        jira_ticket, library_id,
    )

    results_file_resource, results_file_instance = tantalus_api.add_file(
        results_storage_name, results_filepath, update=rerun)
    results_file_pk = results_file_resource['id']

    metadata_file_resource, metadata_file_instance = tantalus_api.add_file(
        results_storage_name, metadata_filepath, update=rerun)
    metadata_file_pk = metadata_file_resource['id']

    results = tantalus_api.get_or_create(
        'results',
        name=results_name,
        results_type=results_type,
        results_version=results_version,
        libraries=[library_pk],
        analysis=analysis['id'],
        file_resources=[results_file_pk, metadata_file_pk],
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
@click.argument('results_storage_name', nargs=1)
@click.option('--archive_storage_name', required=False)
@click.option('--jira_ticket', multiple=True)
@click.option('--rerun', is_flag=True)
def run_all_analyses(results_storage_name, archive_storage_name=None, jira_ticket=None, rerun=False):
    tantalus_api = dbclients.tantalus.TantalusApi()

    datasets = get_unprocessed_hmmcopy(tantalus_api, hmmcopy_tickets=jira_ticket, rerun=rerun)

    logging.info('processing {} datasets'.format(len(datasets)))

    for jira_ticket, dataset in datasets.items():
        logging.info('processing dataset {}, ticket {}'.format(dataset['name'], jira_ticket))

        try:
            run_analysis(
                tantalus_api, dataset, jira_ticket, results_storage_name,
                archive_storage_name=archive_storage_name, rerun=rerun)
        except Exception as e:
            logging.exception('processing of dataset {} failed with exception {}'.format(dataset['name'], e))


if __name__ == "__main__":
    # Set up the root logger
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    run_all_analyses()

