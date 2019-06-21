import collections
import logging
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from dbclients.basicclient import NotFoundError


def search_hmmcopy_analyses(
        tantalus_api,
        library_id,
        aligner_name='BWA_ALN_0_5_7',
        retain_incomplete=False,
        require_unique=True,
    ):
    """ Search for hmmcopy analyses for a specific library.
    """
    logging.info('hmmcopy data for {}'.format(library_id))

    analyses = list(tantalus_api.list(
        'analysis',
        analysis_type__name='hmmcopy',
        input_datasets__library__library_id=library_id,
    ))

    filtered_analyses = []
    for analysis in analyses:
        aligners = set()
        is_complete = True
        for dataset_id in analysis['input_datasets']:
            dataset = tantalus_api.get('sequencedataset', id=dataset_id)
            if not dataset['is_complete']:
                is_complete = False
            aligners.add(dataset['aligner'])
        if not retain_incomplete and not is_complete:
            continue
        if len(aligners) != 1:
            raise Exception('found {} aligners for analysis {}'.format(
                len(aligners), analysis['id']))
        aligner = aligners.pop()
        if aligner == aligner_name:
            filtered_analyses.append(analysis)

    if require_unique and len(filtered_analyses) != 1:
        raise Exception('found {} hmmcopy analyses for {}: {}'.format(
            len(filtered_analyses), library_id, [a['id'] for a in filtered_analyses]))

    if require_unique:
        return filtered_analyses[0]

    return filtered_analyses


def search_cell_cycle_results(
        tantalus_api,
        hmmcopy_ticket_id,
        version='v0.0.2',
        results_storage_name='singlecellblob_results',
):
    """ Import cell cycle predictions for a list of libraries
    """
    hmmcopy_analysis = tantalus_api.get(
        'analysis',
        analysis_type__name='hmmcopy',
        jira_ticket=hmmcopy_ticket_id,
    )

    hmmcopy_results = tantalus_api.get(
        'resultsdataset',
        analysis=hmmcopy_analysis['id'],
    )

    assert len(hmmcopy_results['libraries']) == 1
    library_id = hmmcopy_results['libraries'][0]['library_id']

    logging.info('cell cycle data for {}'.format(library_id))

    analysis = tantalus_api.get(
        'analysis',
        analysis_type='cell_state_classifier',
        version=version,
        input_results__id=hmmcopy_results['id'],
    )

    results = tantalus_api.get(
        'resultsdataset',
        analysis=analysis['id'],
    )

    return results


def search_image_feature_results(
        tantalus_api,
        library_id,
        version='v0.0.1',
        results_storage_name='singlecellblob_results',
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

