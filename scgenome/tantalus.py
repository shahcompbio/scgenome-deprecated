import sys
import logging
import dbclients
import scgenome.utils
from scgenome.constants import LOGGING_FORMAT, CACHE_DIR, PROPS_LENGTH, \
    CELL_ID, ORIGIN_ID
import os
import numpy as np
import pandas as pd

from scgenome.loaders.qc import load_cached_qc_data
from scgenome.db.qc import cache_qc_results
from scgenome.analyses.infer_clones import load_cell_cycle_data
from scgenome.qc import qc_cn


def get_local_cache_dir():
    home = os.path.expanduser("~")
    cache = os.path.join(home, CACHE_DIR)
    if not os.path.exists(cache):
        os.makedirs(cache)
    return cache


def get_data(hmmcopy_tickets, sample_ids, cached=False,
             local_cache_directory=get_local_cache_dir()):
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

    tantalus_api = dbclients.tantalus.TantalusApi()

    cn_data = []
    segs_data = []
    metrics_data = []
    align_metrics_data = []
    for jira_ticket in hmmcopy_tickets:
        analysis = tantalus_api.get('analysis', analysis_type__name='hmmcopy',
                                    jira_ticket=jira_ticket)

        if not cached:
            cache_qc_results(jira_ticket, local_cache_directory)
        hmmcopy_data = load_cached_qc_data(jira_ticket, local_cache_directory,
                                           sample_ids=sample_ids)

        cn_data.append(hmmcopy_data['hmmcopy_reads'])
        segs_data.append(hmmcopy_data['hmmcopy_segs'])
        metrics_data.append(hmmcopy_data['hmmcopy_metrics'])
        align_metrics_data.append(hmmcopy_data['align_metrics'])

    cn_data = scgenome.utils.concat_with_categories(cn_data)
    segs_data = scgenome.utils.concat_with_categories(segs_data)
    metrics_data = scgenome.utils.concat_with_categories(metrics_data)
    align_metrics_data = scgenome.utils.concat_with_categories(align_metrics_data)

    if 'is_s_phase' not in metrics_data:
        cell_cycle_data = load_cell_cycle_data(
            tantalus_api,
            analysis['jira_ticket'])
        cell_cycle_data['cell_id'] = cell_cycle_data['cell_id'].astype('category')

        scgenome.utils.union_categories([
            cn_data,
            metrics_data,
            align_metrics_data,
            cell_cycle_data])

        metrics_data = metrics_data.merge(cell_cycle_data, how='left')
        assert 'cell_id_x' not in metrics_data

    # TODO return dictionary instead so we know which value were unpacking
    return cn_data, segs_data, metrics_data, align_metrics_data

def spike_in(num_cells, hmmcopy_tickets, sample_ids, cached=False,
             local_cache_directory=get_local_cache_dir(), proportions=None,
             id_field_name=CELL_ID, origin_field_name=ORIGIN_ID, seed=None):
    """

   :param hmmcopy_tickets: list of hmmcopy tickets
   :param sample_ids: list of lists of sample ids. first element corresponds
   to first element in `hmmcopy_tickets` etc.
   :param cached:
   :param local_cache_directory:
   :return:
   """
    if len(hmmcopy_tickets) != len(sample_ids):
        raise ValueError('hmmcopy_tickets and sample_ids must be same length')
    if proportions is not None:
        if len(proportions) != len(hmmcopy_tickets):
            raise ValueError(PROPS_LENGTH)
    else:
        if seed is not None:
            np.random.seed(seed)
        proportions = np.random.normal(size=len(sample_ids))
        proportions = proportions / proportions.sum()

    cell_counts = num_cells * proportions

    datasets = []
    sub_datasets = []
    for i in range(len(hmmcopy_tickets)):
        data = get_data(hmmcopy_tickets[i], sample_ids[i], cached,
                        local_cache_directory)
        cn, cn_data = qc_cn(data[2], data[0])
        cn_data[origin_field_name] = i

        sub_cn_data = subsample_cn_data(cn_data, cell_counts[i], id_field_name,
                                        seed)
        datasets.append(cn_data)
        sub_datasets.append(sub_cn_data)

    mixed = pd.concatenat(sub_datasets)
    return {"mixed_cn_data": mixed, "full_datasets": datasets,
            "proportions": proportions, "cell_counts": cell_counts}


def subsample_cn_data(cn_data, num_cells, id_field_name=CELL_ID,
                      seed=None):
    if seed is not None:
        np.random.seed(seed)

    keep_ids = cn_data[id_field_name].unique().sample(num_cells,
                                                      random_state=seed)

    return cn_data[cn_data[id_field_name] in keep_ids]