import collections
import logging
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from dbclients.basicclient import NotFoundError


def _get_analysis_inputs_info(tantalus_api, analysis):
    aligners = set()
    is_complete = True
    bam_datasets = list(tantalus_api.list('sequencedataset', analysis__jira_ticket=analysis['jira_ticket'], dataset_type='BAM'))
    for dataset in bam_datasets:
        if not dataset['is_complete']:
            is_complete = False
        aligners.add(dataset['aligner'])
    if len(aligners) != 1:
        logging.warning('found {} aligners for analysis {}'.format(
            len(aligners), analysis['id']))
        return None
    aligner = aligners.pop()

    info = {
        'is_complete': is_complete,
        'aligner': aligner,
    }

    return info


def search_hmmcopy_analysis(
        tantalus_api,
        library_id,
        aligner_name='BWA_ALN_0_5_7',
    ):
    """ Search for hmmcopy results and analyses for a specific library.
    """
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
        analysis = tantalus_api.get('analysis', id=results['analysis'])

        results_analysis[results['id']] = analysis

        info = _get_analysis_inputs_info(tantalus_api, analysis)

        if info is None:
            continue

        info['results_id'] = results['id']
        info['jira_ticket'] = analysis['jira_ticket']
        info['jira_ticket_number'] = int(analysis['jira_ticket'].split('-')[1])

        results_info.append(info)

    results_info = pd.DataFrame(results_info)

    if aligner_name is not None:
        results_info = results_info.query(f'aligner == "{aligner_name}"')

    if len(results_info) == 0:
        raise ValueError(f'no results for library {library_id} with aligner {aligner_name}')

    # If any of the results are complete, select amongst those
    if results_info['is_complete'].any():
        results_info = results_info.query('is_complete')

    else:
        logging.warning(f'no complete results for library {library_id}')

    # Select the latest by jira ticket number for multiple
    if len(results_info) > 1:
        logging.warning(f'selecting latest of {len(results_info)} results for {library_id} by jira ticket number')
    results_id = results_info.sort_values('jira_ticket_number').iloc[-1, :]['results_id']

    return results_analysis[results_id]


def search_cell_cycle_results(
        tantalus_api,
        hmmcopy_results,
        version='v0.0.1',
        results_storage_name='singlecellresults',
):
    """ Import cell cycle predictions for a list of libraries
    """
    logging.info(f'cell cycle data for {library_id}')

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

