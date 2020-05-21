import collections
import logging
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from dbclients.basicclient import NotFoundError


def get_analysis_lanes(tantalus_api, analysis):
    bam_datasets = list(tantalus_api.list('sequencedataset', analysis__jira_ticket=analysis['jira_ticket'], dataset_type='BAM'))
    lane_set = set()
    for dataset in bam_datasets:
        for lane in dataset['sequence_lanes']:
            lane_set.add('{}_{}'.format(lane['flowcell_id'], lane['lane_number']))
    return lane_set


def get_analysis_inputs_info(tantalus_api, analysis):
    assert analysis['analysis_type'] in ('qc', 'hmmcopy')

    if analysis['analysis_type'] == 'qc':
        bam_datasets = list(tantalus_api.list(
            'sequencedataset',
            analysis__jira_ticket=analysis['jira_ticket'],
            dataset_type='BAM',
        ))

    elif analysis['analysis_type'] == 'hmmcopy':
        bam_datasets = []
        for dataset_id in analysis['input_datasets']:
            bam_datasets.append(tantalus_api.get('sequencedataset', id=dataset_id))

    if len(bam_datasets) == 0:
        logging.error(f'no datasets for analysis {analysis["id"]}')
        return

    field_names = [
        'is_complete',
        'reference_genome',
        'aligner',
        'lane_set',
    ]
    
    info = {}
    for field_name in field_names:
        field_values = set()
        for dataset in bam_datasets:
            if field_name == 'lane_set':
                lanes = []
                for lane in dataset['sequence_lanes']:
                    lanes.append('{}_{}'.format(lane['flowcell_id'], lane['lane_number']))
                lanes = sorted(lanes)
                lanes = ','.join(lanes)
                field_values.add(lanes)

            else:
                field_values.add(dataset[field_name])

        if len(field_values) != 1:
            raise ValueError(f'found {field_values} for {field_name} for analysis {analysis["id"]}')
        
        info[field_name] = field_values.pop()

    return info


def search_hmmcopy_analysis(
        library_id,
        aligner_name='BWA_ALN',
        reference_genome='HG19',
    ):
    """ Search for hmmcopy results and analyses for a specific library.
    """
    tantalus_api = dbclients.tantalus.TantalusApi()

    logging.info('searching for hmmcopy data for {}'.format(library_id))

    library_results = list(tantalus_api.list(
        'resultsdataset',
        results_type='hmmcopy',
        libraries__library_id=library_id,
    ))

    if len(library_results) == 0:
        raise ValueError(f'no results for library {library_id}')

    results_analysis = {}

    results_info = []

    for results in library_results:
        logging.info('found hmmcopy data {}'.format(results['name']))

        analysis = tantalus_api.get('analysis', id=results['analysis'])

        results_analysis[results['id']] = analysis

        info = get_analysis_inputs_info(tantalus_api, analysis)

        if info is None:
            continue

        info['results_id'] = results['id']
        info['jira_ticket'] = analysis['jira_ticket']
        info['jira_ticket_number'] = int(analysis['jira_ticket'].split('-')[1])

        results_info.append(info)

    results_info = pd.DataFrame(results_info)

    if aligner_name is not None:
        results_info = results_info[results_info['aligner'].str.startswith(aligner_name)]

    if len(results_info) == 0:
        raise ValueError(f'no results for library {library_id} with aligner {aligner_name}')

    if reference_genome is not None:
        results_info = results_info.query(f'reference_genome == "{reference_genome}"')

    if len(results_info) == 0:
        raise ValueError(f'no results for library {library_id} with reference_genome {reference_genome}')

    # If any of the results are complete, select amongst those
    if results_info['is_complete'].any():
        results_info = results_info.query('is_complete')

    else:
        logging.warning(f'no complete results for library {library_id}')

    # Select the most recent by jira ticket number for multiple
    if len(results_info) > 1:
        logging.warning(f'selecting most recent of {len(results_info)} results for {library_id} by jira ticket number')
    results_info = results_info.sort_values('jira_ticket_number', ascending=False)

    results_id = results_info.iloc[0, :]['results_id']
    jira_ticket = results_info.iloc[0, :]['jira_ticket']

    logging.info(f'found results from jira ticket {jira_ticket}') 

    return results_analysis[results_id]


def search_cell_cycle_results(
        tantalus_api,
        hmmcopy_results,
        version='v0.0.1',
        results_storage_name='singlecellresults',
):
    """ Import cell cycle predictions for a list of libraries
    """
    analysis = tantalus_api.get(
        'analysis',
        analysis_type='cell_state_classifier',
        version=version,
        input_results__id=hmmcopy_results['id'],
    )

    results = tantalus_api.get(
        'resultsdataset',
        analysis__jira_ticket=analysis['id'],
    )

    return results, analysis


def search_image_feature_results(
        tantalus_api,
        library_id,
        version='v0.0.1',
        results_storage_name='singlecellresults',
):
    """ Import image features for a list of libraries
    """
    features = tantalus_api.get(
        'results',
        results_type='CELLENONE_FEATURES',
        results_version=version,
        libraries__library_id=library_id,
    )

    return features


def search_image_results(
        tantalus_api,
        library_id,
        version='v1',
        results_storage_name='singlecellresults',
):
    """ Import image features for a list of libraries
    """
    features = tantalus_api.get(
        'results',
        results_type='CELLENONE_IMAGES',
        results_version=version,
        libraries__library_id=library_id,
    )

    return features

