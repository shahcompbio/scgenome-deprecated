import collections
import logging
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from dbclients.basicclient import NotFoundError


def search_hmmcopy_results(
        tantalus_api,
        library_id,
        aligner_name='BWA_ALN_0_5_7',
):
    """ Search for hmmcopy results and analyses for a specific library.
    """
    logging.info('hmmcopy data for {}'.format(library_id))

    analyses = list(tantalus_api.list(
        'analysis',
        analysis_type__name='hmmcopy',
        input_datasets__library__library_id=library_id,
    ))
    aln_analyses = []
    for analysis in analyses:
        aligners = set()
        is_complete = True
        for dataset_id in analysis['input_datasets']:
            dataset = tantalus_api.get('sequencedataset', id=dataset_id)
            if not dataset['is_complete']:
                is_complete = False
            aligners.add(dataset['aligner'])
        if len(aligners) != 1:
            # HACK: should be an exception but will remove when datasets are cleaned up
            continue
            raise Exception('found {} aligners for analysis {}'.format(
                len(aligners), analysis['id']))
        aligner = aligners.pop()
        if is_complete and aligner == aligner_name:
            aln_analyses.append(analysis)
    if len(aln_analyses) != 1:
        raise Exception('found {} hmmcopy analyses for {}: {}'.format(
            len(aln_analyses), library_id, [a['id'] for a in aln_analyses]))
    analysis = aln_analyses[0]
    results = tantalus_api.get(
        'resultsdataset',
        analysis=analysis['id'],
    )

    return results, analysis


def search_cell_cycle_results(
        tantalus_api,
        library_id,
        hmmcopy_results,
        version='v0.0.1',
        results_storage_name='singlecellblob_results',
):
    """ Import cell cycle predictions for a list of libraries
    """
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

    return results, analysis


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

