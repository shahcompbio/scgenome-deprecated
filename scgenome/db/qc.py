
import dbclients.tantalus
import datamanagement.transfer_files

import scgenome.loaders.align
import scgenome.loaders.hmmcopy


def cache_qc_results(
        ticket_id,
        local_cache_directory,
        full_dataset=False,
        results_storage_name='singlecellresults',
    ):
    tantalus_api = dbclients.tantalus.TantalusApi()

    ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

    results_ids = set()

    for results in ticket_results:
        if full_dataset:
            csv_suffixes = (None,)
        
        else:
            if results['results_type'] == 'align':
                csv_suffixes = scgenome.loaders.align.table_suffixes[results['results_version']]

            elif results['results_type'] == 'hmmcopy':
                csv_suffixes = scgenome.loaders.hmmcopy.table_suffixes[results['results_version']]

            elif results['results_type'] == 'annotation':
                csv_suffixes = scgenome.loaders.annotation.table_suffixes[results['results_version']]

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

        results_ids.add(results['id'])

    for results_id in results_ids:
        filepaths = datamanagement.transfer_files.cache_dataset(
            tantalus_api,
            results_id,
            'resultsdataset',
            results_storage_name,
            local_cache_directory,
            suffix_filter='metadata.yaml',
        )

        if len(filepaths) != 1:
            raise Exception(f'found {len(filepaths)} filepaths for metadata.yaml, results {results["id"]}')


