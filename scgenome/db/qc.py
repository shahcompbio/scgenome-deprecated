import os
import logging
import dbclients.tantalus
import datamanagement.transfer_files

import scgenome.loaders.align
import scgenome.loaders.hmmcopy
import scgenome.loaders.qc


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
            csv_suffixes = ((None,None),)
        
        else:
            if results['results_type'] == 'alignment':
                csv_suffixes = scgenome.loaders.align.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'hmmcopy':
                csv_suffixes = scgenome.loaders.hmmcopy.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'annotation':
                csv_suffixes = scgenome.loaders.annotation.table_suffixes[results['results_version']] + ((None, 'metadata.yaml'),)

            elif results['results_type'] == 'cell_state_prediction':
                csv_suffixes = ((None,None),)

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
        do_caching=False,
    ):

    results_tables = {}

    for ticket_id in ticket_ids:
        if do_caching:
            cache_qc_results(ticket_id, local_directory)

        ticket_directory = os.path.join(local_directory, ticket_id)
        ticket_results = scgenome.loaders.qc.load_qc_data(ticket_directory, sample_ids=sample_ids)

        for table_name, table_data in ticket_results.items():
            if table_name not in results_tables:
                results_tables[table_name] = []
            results_tables[table_name].append(table_data)

    for table_name, table_data in results_tables.items():
        results_tables[table_name] = scgenome.utils.concat_with_categories(table_data)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables
