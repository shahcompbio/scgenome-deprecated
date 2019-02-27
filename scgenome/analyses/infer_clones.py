import sys
import os
import logging
import click
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wget
import functools
import itertools
import pickle

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
        sample_ids=sample_ids,
        #ploidy_solution='2',
        #subsample=0.25,
    )

    cn_data = results['hmmcopy_reads']
    metrics_data = results['hmmcopy_metrics']

    cell_cycle_data = get_cell_cycle_data(tantalus_api, library_ids, hmmcopy_results)
    cell_cycle_data['cell_id'] = pd.Categorical(cell_cycle_data['cell_id'], categories=metrics_data['cell_id'].cat.categories)
    metrics_data = metrics_data.merge(cell_cycle_data)

    image_feature_data = get_image_feature_data(tantalus_api, library_ids)

    # Read count filtering
    metrics_data = metrics_data[metrics_data['total_mapped_reads'] > 500000]

    # Filter by experimental condition
    metrics_data = metrics_data[~metrics_data['experimental_condition'].isin(['NTC'])]

    cell_ids = metrics_data['cell_id']
    cn_data = cn_data[cn_data['cell_id'].isin(cell_ids)]

    return cn_data, metrics_data, image_feature_data


def infer_clones(cn_data, metrics_data, results_prefix):
    total_cells = len(metrics_data.index)

    # Remove s phase cells
    metrics_data = metrics_data[~metrics_data['is_s_phase']]

    # Remove low quality cells
    metrics_data = metrics_data[metrics_data['quality'] > 0.5]

    # Remove low coverage cells
    metrics_data = metrics_data[metrics_data['total_mapped_reads'] > 500000]

    logging.info('filtering {} of {} cells'.format(
        len(metrics_data.index), total_cells))
    cn_data = cn_data.merge(metrics_data[['cell_id']].drop_duplicates())
    assert isinstance(cn_data['cell_id'].dtype, pd.api.types.CategoricalDtype)

    logging.info('creating copy number matrix')
    cn = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level='cell_id').fillna(0)
    )

    logging.info('clustering copy number')
    clusters = scgenome.cncluster.umap_hdbscan_cluster(cn)

    logging.info('merging clusters')
    cn_data = cn_data.merge(clusters[['cell_id', 'cluster_id']].drop_duplicates())

    logging.info('plotting clusters to {}*'.format(results_prefix + '_initial'))
    plot_clones(cn_data, 'cluster_id', results_prefix + '_initial')

    clone_cn_data = (
        cn_data
            .groupby(['chr', 'start', 'end', 'cluster_id'])
            .agg({'copy': np.mean, 'state': np.median})
            .reset_index()
    )
    clone_cn_data['state'] = clone_cn_data['state'].round().astype(int)

    return clusters, clone_cn_data


def reassign_cells(cn_data, clone_cn_data, results_prefix):
    """
    """
    logging.info('Create matrix of cn data for all cells')
    cell_cn_matrix = (
        cn_data
            .set_index(['chr', 'start', 'end', 'cell_id'])['copy']
            .unstack(level=['cell_id']).fillna(0)
    )

    logging.info('Create a matrix of cn data for filtered clones')
    clone_cn_matrix = (
        clone_cn_data
            .set_index(['chr', 'start', 'end', 'cluster_id'])['state']
            .unstack(level=['cluster_id']).fillna(0)
    )

    logging.info('Calculating clone cell correlation')
    cell_clone_corr = {}
    for cluster_id in clone_cn_matrix.columns:
        logging.info('Calculating correlation for clone {}'.format(cluster_id))
        cell_clone_corr[cluster_id] = cell_cn_matrix.corrwith(clone_cn_matrix[cluster_id])

    reclusters = pd.DataFrame(cell_clone_corr).idxmax(axis=1).dropna().astype(int)
    reclusters.name = 'recluster_id'
    reclusters = reclusters.reset_index()

    logging.info('merging clusters')
    cn_data = cn_data.merge(reclusters[['cell_id', 'recluster_id']].drop_duplicates())

    logging.info('plotting clusters to {}*'.format(results_prefix + '_recluster'))
    plot_clones(cn_data, 'recluster_id', results_prefix + '_recluster')

    return [reclusters]


def plot_clones(cn_data, cluster_col, plots_prefix):
    plot_data = cn_data.copy()
    bin_filter = (plot_data['gc'] <= 0) | (plot_data['copy'].isnull())
    plot_data.loc[bin_filter, 'state'] = 0
    plot_data.loc[plot_data['copy'] > 5, 'copy'] = 5.
    plot_data.loc[plot_data['copy'] < 0, 'copy'] = 0.

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'copy', cluster_field_name=cluster_col, raw=True)
    fig.savefig(plots_prefix + '_raw_cn.pdf')

    fig = plt.figure(figsize=(20, 30))
    matrix_data = scgenome.cnplot.plot_clustered_cell_cn_matrix_figure(
        fig, plot_data, 'state', cluster_field_name=cluster_col)
    fig.savefig(plots_prefix + '_cn_state.pdf')


# Until we switch to python3 and can use pickle, cache using hdf5
# and assume a list of pandas dataframes
def memoize(func):
    @functools.wraps(func)
    def memoized_func(cache_filename, *args, **kwargs):
        if os.path.exists(cache_filename):
            logging.info('reading existing data from {}'.format(cache_filename))
            with open(cache_filename) as f:
                data = []
                with pd.HDFStore(cache_filename) as store:
                    for i in itertools.count():
                        key = 'table_{}'.format(i)
                        if key not in store:
                            break
                        data.append(store.get('table_{}'.format(i)))
        else:
            data = func(*args, **kwargs)
            logging.info('writing data to {}'.format(cache_filename))
            with pd.HDFStore(cache_filename, 'w') as store:
                for i, table in enumerate(data):
                    store.put('table_{}'.format(i), table, format='table')
        return data
    return memoized_func


retrieve_data_with_cache = memoize(retrieve_data)
infer_clones_with_cache = memoize(infer_clones)
reassign_cells_with_cache = memoize(reassign_cells)


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

    logging.info('retrieving cn data')
    cn_data, metrics_data, image_feature_data = retrieve_data_with_cache(
        raw_data_filename,
        tantalus_api,
        library_ids,
        sample_ids,
        local_storage_directory,
        results_prefix,
    )

    clones_filename = results_prefix + '_clones.h5'
    logging.info('inferring clones')
    shape_check = cn_data.shape
    logging.info('cn_data shape {}'.format(shape_check))
    clusters, clone_cn_data = infer_clones_with_cache(
        clones_filename,
        cn_data,
        metrics_data,
        results_prefix,
    )
    assert cn_data.shape == shape_check

    reassign_clones_filename = results_prefix + '_reassign_clones.h5'
    reclusters = reassign_cells_with_cache(
        reassign_clones_filename,
        cn_data,
        clone_cn_data,
        results_prefix,
    )

    print(cn_data.head())
    print(metrics_data.head())
    print(image_feature_data.head())
    print(clone_cn_data.head())
    print(clusters.head())


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    infer_clones_cmd()
