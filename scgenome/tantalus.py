import sys
import logging
import dbclients
import scgenome.utils
from cnv_cluster.constants import LOGGING_FORMAT, CACHE_DIR
import os

from scgenome.loaders.qc import load_cached_qc_data
from scgenome.db.qc import cache_qc_results
from scgenome.analyses.infer_clones import load_cell_cycle_data


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

    return cn_data, segs_data, metrics_data, align_metrics_data