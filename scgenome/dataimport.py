import collections
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files


default_hmmcopy_reads_cols = (
    'chr',
    'start',
    'end',
    'cell_id',
    'gc',
    'reads',
    'copy',
    'state',
    'integer_copy_number',
)


def import_cn_data(
        tickets,
        sample_ids,
        local_storage_name,
        cols=default_hmmcopy_reads_cols,
        local_storage_directory=None,
        results_storage_name='singlecellblob_results',
):
    tantalus_api = dbclients.tantalus.TantalusApi()

    results_storage = tantalus_api.get("storage", name=results_storage_name)
    local_storage = tantalus_api.get("storage", name=local_storage_name)

    local_results_client = tantalus_api.get_storage_client(local_storage_name)

    if local_storage_directory is not None:
        tantalus_api.get_or_create(
            "storage_server",
            storage_type="server",
            name=local_storage_name,
            storage_directory=local_storage_directory,
        )

    else:
        tantalus_api.get(
            "storage_server",
            name=local_storage_name,
        )

    results_tables = collections.defaultdict(list)

    results_info = [
        ('{}_hmmcopy', '_hmmcopy.h5', [
            ('hmmcopy_reads', '/hmmcopy/reads/0', default_hmmcopy_reads_cols),
            ('hmmcopy_metrics', '/hmmcopy/metrics/0', None),
        ]),
        ('{}_align', '_alignment_metrics.h5', [
            ('align_metrics', '/alignment/metrics', None),
        ]),
    ]

    for analysis_template, h5_suffix, table_info in results_info:

        for ticket in tickets:
            results = tantalus_api.get('results', name=analysis_template.format(ticket))
    
            datamanagement.transfer_files.transfer_results_dataset(
                tantalus_api,
                results['id'],
                results_storage,
                local_storage,
            )
    
            for file_resource_id in results['file_resources']:
                file_resource = tantalus_api.get('file_resource', id=file_resource_id)

                if file_resource['filename'].endswith(h5_suffix):
                    for table_name, table_key, columns in table_info:
                        filepath = local_results_client.get_url(file_resource['filename'])
                        results_tables[table_name].append(pd.read_hdf(filepath, table_key, columns=columns))

    for table_name, results_data in results_tables.iteritems():
        results_data = pd.concat(results_data, ignore_index=True)
        if 'sample_id' not in results_data:
            results_data['sample_id'] = results_data['cell_id'].apply(lambda a: a.split('-')[0])
        if 'library_id' not in results_data:
            results_data['library_id'] = results_data['cell_id'].apply(lambda a: a.split('-')[1])
        results_data = results_data[results_data['sample_id'].isin(sample_ids)]
        results_tables[table_name] = results_data

    return results_tables


