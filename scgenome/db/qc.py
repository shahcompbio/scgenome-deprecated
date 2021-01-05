import logging
import os

import datamanagement.transfer_files
import dbclients.tantalus
import scgenome.loaders.align
import scgenome.loaders.hmmcopy
import scgenome.loaders.qc
from scgenome.db.qc_from_files import _aggregate_results_tables
from scgenome.db.qc_from_files import _concat_results_tables


def cache_qc_results(
        ticket_id,
        local_cache_directory,
        full_dataset=False,
        results_storage_name='singlecellresults',
    ):
    tantalus_api = dbclients.tantalus.TantalusApi()

    ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

    for results in ticket_results:
        logging.info(f'found results {results["id"]} with type {results["results_type"]} for ticket {ticket_id}')

        if full_dataset:
            csv_suffixes = ((None, None),)

        else:
            if results['results_type'] == 'alignment':
                csv_suffixes = scgenome.loaders.align.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'hmmcopy':
                csv_suffixes = scgenome.loaders.hmmcopy.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'annotation':
                csv_suffixes = scgenome.loaders.annotation.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'cell_state_prediction':
                csv_suffixes = ((None, None),)

            else:
                continue

        for _, csv_suffix in csv_suffixes:
            filepaths = datamanagement.transfer_files.cache_dataset(
                tantalus_api,
                results['id'],
                'resultsdataset',
                results_storage_name,
                local_cache_directory,
                suffix_filter=csv_suffix,
            )

            if csv_suffix is not None and len(filepaths) != 1:
                raise Exception(f'found {len(filepaths)} filepaths for {csv_suffix}, results {results["id"]}')


def get_qc_data(
        ticket_ids,
        local_directory,
        sample_ids=None,
        additional_hmmcopy_reads_cols=None,
        do_caching=False,
    ):
    results_tables = {}

    for ticket_id in ticket_ids:
        if do_caching:
            cache_qc_results(ticket_id, local_directory)

        ticket_directory = os.path.join(local_directory, ticket_id)

        ticket_results = scgenome.loaders.qc.load_qc_data(
            ticket_directory, sample_ids=sample_ids,
            additional_hmmcopy_reads_cols=additional_hmmcopy_reads_cols)

        results_tables = _aggregate_results_tables(results_tables, ticket_results)

    results_tables = _concat_results_tables(results_tables)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables
