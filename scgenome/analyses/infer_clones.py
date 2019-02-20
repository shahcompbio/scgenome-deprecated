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
import scgenome.cncluster
import scgenome.cnplot

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
    cell_cycle_data = pd.concat(cell_cycle_data, ignore_index=True, sort=True)

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
    image_feature_data = pd.concat(image_feature_data, ignore_index=True, sort=True)

    return image_feature_data


def retrieve_data(tantalus_api, library_ids, sample_ids, local_storage_directory, results_prefix):
    hmmcopy_results, hmmcopy_tickets = get_hmmcopy_analyses(tantalus_api, library_ids)

    results = scgenome.dataimport.import_cn_data(
        hmmcopy_tickets,
        local_storage_directory,
        #ploidy_solution='2',
        #subsample=0.25,
    )

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']

    metrics_data['sample_id'] = [a.split('-')[0] for a in metrics_data.index.get_level_values('cell_id')]
    metrics_data['library_id'] = [a.split('-')[1] for a in metrics_data.index.get_level_values('cell_id')]

    cell_cycle_data = get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)
    cell_cycle_data.set_index('cell_id', inplace=True)
    metrics_data = metrics_data.merge(cell_cycle_data, left_index=True, right_index=True)

    image_feature_data = get_image_feature_data(tantalus_api, library_ids)

    # Sample filtering
    metrics_data = metrics_data[metrics_data['sample_id'].isin(sample_ids)]

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads'] > 500000]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data.index.get_level_values('cell_id')
    cn_data = cn_data[cn_data.index.get_level_values('cell_id').isin(cell_ids)]

    return cn_data, metrics_data, image_feature_data


def infer_clones(cn_data, metrics_data):
    total_cells = len(metrics_data.index)

    # Remove s phase cells
    metrics_data = metrics_data[~metrics_data['is_s_phase']]

    # Remove low quality cells
    metrics_data = metrics_data[metrics_data['quality'] > 0.5]

    # Remove low coverage cells
    metrics_data = metrics_data[metrics_data['total_mapped_reads'] > 500000]

    logging.info('filtering {} of {} cells'.format(
        len(metrics_data.index), total_cells))

    cn_data = cn_data[cn_data.index.get_level_values('cell_id').isin(
        metrics_data.index.get_level_values('cell_id'))]

    logging.info('creating copy number matrix')

    cn = cn_data['copy'].unstack(level='cell_id').fillna(0)

    print cn.head()

    cluster_df = scgenome.cncluster.umap_hdbscan_cluster(cn)

    print cluster_df.head()

    cn_data.set_index(
        cluster_df.set_index('cell_id').loc[
            cn_data.index.get_level_values('cell_id'), 'cluster_id'],
        append=True, inplace=True)
    clone_cn = (
        cn_data.groupby(level=['chr', 'start', 'end', 'cluster_id'])['copy']
        .median().rename('clone_cn').reset_index())

    print cluster_df.groupby('cluster_id').size().rename('size').reset_index()

    return cluster_df, clone_cn


def plot_clones(results_prefix, cn_data):
    plot_data = cn_data.copy()
    bin_filter = (plot_data['gc'] <= 0) | (plot_data['copy'].isnull())
    plot_data.loc[bin_filter, 'state'] = 0
    plot_data.loc[plot_data['copy'] > 5, 'copy'] = 5.
    plot_data.loc[plot_data['copy'] < 0, 'copy'] = 0.

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'copy', raw=True)
    fig.savefig(results_prefix + '_raw_cn.pdf')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'state')
    fig.savefig(results_prefix + '_cn_state.pdf')





@click.command()
@click.argument('library_ids_filename')
@click.argument('sample_ids_filename')
@click.argument('results_prefix')
@click.argument('local_storage_directory')
def infer_clones_cmd(library_ids_filename, sample_ids_filename, results_prefix, local_storage_directory):
    tantalus_api = dbclients.tantalus.TantalusApi()

    library_ids = [l.strip() for l in open(library_ids_filename).readlines()]
    sample_ids = [l.strip() for l in open(sample_ids_filename).readlines()]

    raw_data_filename = results_prefix + '_raw_data.h5'

    if os.path.exists(raw_data_filename):
        logging.info('reading previous data')

        cn_data = pd.read_hdf(raw_data_filename, 'cn_data')
        metrics_data = pd.read_hdf(raw_data_filename, 'metrics')
        image_feature_data = pd.read_hdf(raw_data_filename, 'image_features')

    else:
        logging.info('retrieving data')

        cn_data, metrics_data, image_feature_data = retrieve_data(
            tantalus_api, library_ids, sample_ids, local_storage_directory, results_prefix)

        cn_data.to_hdf(raw_data_filename, 'cn_data')
        metrics_data.to_hdf(raw_data_filename, 'metrics')
        image_feature_data.to_hdf(raw_data_filename, 'image_features')

    clones_filename = results_prefix + '_clones.h5'

    if os.path.exists(clones_filename):
        clone_cn = pd.read_hdf(clones_filename, 'clone_cn')
        clusters = pd.read_hdf(clones_filename, 'clusters')

    else:
        logging.info('inferring clones')

        shape_check = cn_data.shape
        logging.info('cn_data shape {}'.format(shape_check))
        clusters, clone_cn = infer_clones(cn_data, metrics_data)
        assert cn_data.shape == shape_check

        clone_cn.to_hdf(clones_filename, 'clone_cn')
        clusters.to_hdf(clones_filename, 'clusters')

    plot_clones(results_prefix, cn_data)

    print cn_data.head()
    print metrics_data.head()
    print image_feature_data.head()
    print clone_cn.head()
    print clusters.head()


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones_cmd()
