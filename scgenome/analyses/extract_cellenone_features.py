import os
import sys
import click
import logging
import itertools
import pandas as pd

import dbclients.tantalus
import dbclients.colossus
import datamanagement.transfer_files
from datamanagement.utils.utils import make_dirs

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"



analysis_type = 'cellenone_extraction'
analysis_version = 'v0.0.1'


results_type = 'CELLENONE_FEATURES'
results_version = 'v0.0.1'


txt_cols = [
    'YPos',
    'XPos',
    'Circ', 
    'Elong',
    'Diameter',
]


def calculate_brief_cell_id(data, swap):
    if swap:
        data['brief_cell_id'] = data.apply(
            lambda row: '{}-R{:02d}-C{:02d}'.format(
                row['library_id'], int(row['XPos']), int(row['YPos'])), axis=1)

    else:
        data['brief_cell_id'] = data.apply(
            lambda row: '{}-R{:02d}-C{:02d}'.format(
                row['library_id'], int(row['YPos']), int(row['XPos'])), axis=1)


def calculate_swap(data, cell_ids, library_id):
    swap_scores = []
    for swap in (True, False):
        calculate_brief_cell_id(data, swap)
        score = data['brief_cell_id'].isin(cell_ids).mean()
        swap_scores.append((score, swap))

        logging.info('library: {}, swap: {}, score: {}'.format(
            library_id, swap, score))

    selected_swap = sorted(swap_scores)[-1][1]

    logging.info('library: {}, selected swap: {}'.format(
        library_id, selected_swap))
    
    calculate_brief_cell_id(data, selected_swap)
    data['selected_swap'] = selected_swap

    return data


def process_cellenone_table(tantalus_api, cellenone_dataset, inputs_storage_name):
    storage_client = tantalus_api.get_storage_client(inputs_storage_name)

    assert len(cellenone_dataset['libraries']) == 1
    library_id = cellenone_dataset['libraries'][0]['library_id']

    txt_file_instances = tantalus_api.get_dataset_file_instances(
        cellenone_dataset['id'], 'resultsdataset', inputs_storage_name,
        filters={'filename__endswith': '.txt'})
    xls_file_instances = tantalus_api.get_dataset_file_instances(
        cellenone_dataset['id'], 'resultsdataset', inputs_storage_name,
        filters={'filename__endswith': '.xls'})

    filenames = set()
    for file_instance in itertools.chain(txt_file_instances, xls_file_instances):
        filename = file_instance['file_resource']['filename']
        if not filename.endswith('.xls') and not filename.endswith('.txt'):
            raise ValueError()
        filenames.add(filename)

    cell_info = dbclients.colossus.get_colossus_sublibraries_from_library_id(library_id, brief=True)

    cell_ids = set(["{}-R{}-C{}".format(library_id, a["row"], a["column"]) for a in cell_info])
    cell_ids = pd.Series(list(cell_ids), name='brief_cell_id')

    cellenone_data = []

    for filename in filenames:
        url = storage_client.get_url(filename)

        if filename.endswith('isolated.xls'):
            data = pd.read_csv(url, sep='\t')
            data['filename_suffix'] = 'isolated.xls'
            data['spotter'] = 'rachael'
        elif filename.endswith('singleprinted.xls'):
            data = pd.read_csv(url, sep='\t')
            data['filename_suffix'] = 'isolated.xls'
            data['spotter'] = 'rachael'
        elif filename.endswith('.txt') and not filename.endswith('Mapping.txt'):
            data = pd.read_csv(url, sep='\t', names=txt_cols)
            data['filename_suffix'] = 'txt'
            data['spotter'] = 'deckard'
        else:
            continue

        # Filter rows for which no cell was found
        # these have misaligned data in the TSV
        data = data[data['XPos'].notnull()]
        data = data[data['YPos'].notnull()]
        data = data[data['Circ'].notnull()]

        if data.empty:
            continue

        data['library_id'] = library_id

        cellenone_data.append(data)

    cellenone_data = pd.concat(cellenone_data, ignore_index=True)

    corrected_data = []
    for spotter, data in cellenone_data.groupby('spotter'):
        if spotter == 'rachael':
            data = calculate_swap(data.copy(), cell_ids, library_id)
        else:
            data = data.copy()
            data['selected_swap'] = False
            calculate_brief_cell_id(data, False)

        corrected_data.append(data)

    cellenone_data = pd.concat(corrected_data, ignore_index=True)

    return cellenone_data


def get_unprocessed_cellenone(tantalus_api):
    datasets = {}

    for results in tantalus_api.list('results', results_type='CELLENONE'):
        assert len(results['libraries']) <= 1

        if len(results['libraries']) == 0:
            logging.warning('no library for {}'.format(results['name']))
            continue
        library_id = results['libraries'][0]['library_id']

        try:
            features = tantalus_api.get(
                'results',
                results_type='CELLENONE_FEATURES',
                results_version=results_version,
                libraries__library_id=library_id,
            )
        except dbclients.basicclient.NotFoundError:
            features = None

        if features is not None:
            logging.info('library {} has features {}'.format(
                library_id, features['name']))
            continue

        assert library_id not in datasets
        datasets[library_id] = results

    return datasets.values()


def run_analysis(
        tantalus_api, jira_ticket, cellenone_dataset,
        inputs_storage_name, results_storage_name,
        archive_storage_name=None):

    results_storage = tantalus_api.get(
        'storage',
        name=results_storage_name,
    )

    assert len(cellenone_dataset['libraries']) == 1
    library_pk = cellenone_dataset['libraries'][0]['id']
    library_id = cellenone_dataset['libraries'][0]['library_id']

    results_filename = os.path.join(
        'single_cell_indexing',
        'Cellenone',
        'Cellenone_data',
        'feature_tables',
        '{results_version}',
        '{library_id}.csv',
    ).format(
        results_version=results_version,
        library_id=library_id,
    )

    results_filepath = os.path.join(
        results_storage['storage_directory'],
        results_filename,
    )

    make_dirs(os.path.dirname(results_filepath))

    cellenone_data = process_cellenone_table(tantalus_api, cellenone_dataset, inputs_storage_name)
    cellenone_data.to_csv(results_filepath)

    analysis_name = '{}_{}_{}'.format(
        analysis_type, analysis_version, library_id)

    analysis = tantalus_api.get_or_create(
        'analysis',
        name=analysis_name,
        analysis_type=analysis_type,
        version=analysis_version,
        jira_ticket=jira_ticket,
        status='complete',
        args={},
        input_results=[cellenone_dataset['id']],
    )

    results_name = '{}_{}_{}'.format(
        results_type, results_version, library_id
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
@click.argument('inputs_storage_name', nargs=1)
@click.argument('results_storage_name', nargs=1)
@click.option('--archive_storage_name', required=False)
def run_all_analyses(jira_ticket, inputs_storage_name, results_storage_name, archive_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    datasets = get_unprocessed_cellenone(tantalus_api)
    for dataset in datasets:
        logging.info('processing dataset {}'.format(dataset['name']))

        try:
            run_analysis(
                tantalus_api, jira_ticket, dataset, inputs_storage_name, results_storage_name,
                archive_storage_name=archive_storage_name)
        except Exception as e:
            logging.warning('processing of dataset {} failed with exception {}'.format(dataset['name'], e))


if __name__ == "__main__":
    # Set up the root logger
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    run_all_analyses()
