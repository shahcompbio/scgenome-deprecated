import sys
import os
import logging
import click
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wget

import seaborn
import numpy as np
import pandas as pd
import pylab
import sklearn.preprocessing

import scgenome
import scgenome.dataimport

import dbclients.tantalus
from dbclients.basicclient import NotFoundError
import datamanagement.transfer_files


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


results_storage_name = 'singlecellblob_results'


def get_hmmcopy_analyses(tantalus_api, library_ids):
    hmmcopy_results = {}
    hmmcopy_tickets = []

    for library_id in library_ids:
        logging.info('hmmcopy data for {}'.format(library_id))

        analyses = list(tantalus_api.list(
            'analysis',
            analysis_type__name='hmmcopy',
            input_datasets__library__library_id=library_id,
        ))
        aln_analyses = []
        for analysis in analyses:
            aligners = set()
            for dataset_id in analysis['input_datasets']:
                dataset = tantalus_api.get('sequencedataset', id=dataset_id)
                aligners.add(dataset['aligner'])
            if len(aligners) != 1:
                raise Exception('found {} aligners for analysis {}'.format(
                    len(aligners), analysis['id']))
            aligner = aligners.pop()
            if aligner == 'BWA_ALN_0_5_7':
                aln_analyses.append(analysis)
        if len(aln_analyses) != 1:
            raise Exception('found {} hmmcopy analyses for {}: {}'.format(
                len(aln_analyses), library_id, [a['id'] for a in aln_analyses]))
        analysis = aln_analyses[0]
        results = tantalus_api.get(
            'resultsdataset',
            analysis=analysis['id'],
        )
        hmmcopy_results[library_id] = results
        hmmcopy_tickets.append(analysis['jira_ticket'])

    return hmmcopy_results, hmmcopy_tickets


def get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results):
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    cell_cycle_data = []

    for library_id in library_ids:
        logging.info('cell cycle data for {}'.format(library_id))

        classifier = tantalus_api.get(
            'analysis',
            analysis_type='cell_state_classifier',
            version='v0.0.1',
            input_results__id=hmmcopy_results[library_id]['id'],
        )
        features = tantalus_api.get(
            'resultsdataset',
            analysis=classifier['id'],
        )
        file_instances = tantalus_api.get_dataset_file_instances(
            features['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f)
            data['library_id'] = library_id
            cell_cycle_data.append(data)
    cell_cycle_data = pd.concat(cell_cycle_data, ignore_index=True)

    return cell_cycle_data


def get_image_feature_data(tantalus_api, library_ids):
    storage_client = tantalus_api.get_storage_client(results_storage_name)

    image_feature_data = []

    for library_id in library_ids:
        try:
            features = tantalus_api.get(
                'results',
                results_type='CELLENONE_FEATURES',
                results_version='v0.0.1',
                libraries__library_id=library_id,
            )
        except NotFoundError:
            logging.info('no image data for {}'.format(library_id))
            continue
        file_instances = tantalus_api.get_dataset_file_instances(
            features['id'], 'resultsdataset', results_storage_name)
        for file_instance in file_instances:
            f = storage_client.open_file(file_instance['file_resource']['filename'])
            data = pd.read_csv(f, index_col=0)
            data['library_id'] = library_id
            image_feature_data.append(data)
    image_feature_data = pd.concat(image_feature_data, ignore_index=True)

    return image_feature_data


def retrieve_data(tantalus_api, library_ids, local_storage_directory, results_prefix):
    hmmcopy_results, hmmcopy_tickets = get_hmmcopy_analyses(tantalus_api, library_ids)

    results = scgenome.dataimport.import_cn_data(
        hmmcopy_tickets,
        local_storage_directory,
        ploidy_solution='2',
        subsample=0.25,
    )

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']

    cell_cycle_data = get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)

    image_feature_data = get_image_feature_data(tantalus_api, library_ids)

    metrics_data = metrics_data.merge(cell_cycle_data)
    cn_data = cn_data.merge(metrics_data[['cell_id', 'is_s_phase']].drop_duplicates())

    return cn_data, metrics_data, image_feature_data


@click.command()
@click.argument('library_ids_filename')
@click.argument('sample_ids_filename')
@click.argument('results_prefix')
@click.argument('local_storage_directory')
def infer_clones(library_ids_filename, sample_ids_filename, results_prefix, local_storage_directory):
    tantalus_api = dbclients.tantalus.TantalusApi()

    library_ids = [l.strip() for l in open(library_ids_filename).readlines()]

    cn_data_filename = results_prefix + '_cn_data.json.gz'
    metrics_data_filename = results_prefix + '_metrics_data.json.gz'
    image_data_filename = results_prefix + '_image_data.json.gz'

    if os.path.exists(cn_data_filename):
        cn_data = pd.read_json(cn_data_filename)
        metrics_data = pd.read_json(metrics_data_filename)
        image_feature_data = pd.read_json(image_data_filename)

    else:
        cn_data, metrics_data, image_feature_data = retrieve_data(
            tantalus_api, library_ids, local_storage_directory, results_prefix)

        cn_data.to_json(cn_data_filename)
        metrics_data.to_json(metrics_data_filename)
        image_feature_data.to_json(image_data_filename)

    print cn_data.head()
    print metrics_data.head()
    print image_feature_data.head()


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones()
