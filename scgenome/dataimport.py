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
        local_cache_directory,
        cols=default_hmmcopy_reads_cols,
        results_storage_name='singlecellblob_results',
):
    tantalus_api = dbclients.tantalus.TantalusApi()

    local_results_client = tantalus_api.get_cache_client(local_cache_directory)

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
    
            datamanagement.transfer_files.cache_dataset(
                tantalus_api,
                results['id'],
                'resultsdataset',
                results_storage_name,
                local_cache_directory,
            )
    
            for file_resource_id in results['file_resources']:
                file_resource = tantalus_api.get('file_resource', id=file_resource_id)

                if file_resource['filename'].endswith(h5_suffix):
                    for table_name, table_key, columns in table_info:
                        filepath = local_results_client.get_url(file_resource['filename'])
                        data = pd.read_hdf(filepath, table_key, columns=columns)
                        results_tables[table_name].append(data)

    for table_name, results_data in results_tables.iteritems():
        results_data = pd.concat(results_data, ignore_index=True)
        results_data['sample_id'] = results_data['cell_id'].apply(lambda a: a.split('-')[0])
        results_data['library_id'] = results_data['cell_id'].apply(lambda a: a.split('-')[1])
        results_data = results_data[results_data['sample_id'].isin(sample_ids)]
        results_tables[table_name] = results_data

    return results_tables


