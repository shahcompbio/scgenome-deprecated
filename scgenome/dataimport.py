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
        local_cache_directory,
        ploidy_solution='0',
        cols=default_hmmcopy_reads_cols,
        results_storage_name='singlecellblob_results',
        subsample=None,
):
    tantalus_api = dbclients.tantalus.TantalusApi()

    local_results_client = tantalus_api.get_cache_client(local_cache_directory)

    results_tables = collections.defaultdict(dict)

    results_info = [
        ('{}_hmmcopy', '_hmmcopy.h5', [
            ('hmmcopy_reads', '/hmmcopy/reads/'+ploidy_solution, default_hmmcopy_reads_cols, ['chr', 'start', 'end', 'cell_id']),
            ('hmmcopy_metrics', '/hmmcopy/metrics/'+ploidy_solution, None, ['cell_id']),
        ]),
        ('{}_align', '_alignment_metrics.h5', [
            ('align_metrics', '/alignment/metrics', None, ['cell_id']),
        ]),
    ]

    for ticket in tickets:
        for analysis_template, h5_suffix, table_info in results_info:
            results = tantalus_api.get('results', name=analysis_template.format(ticket))

            datamanagement.transfer_files.cache_dataset(
                tantalus_api,
                results['id'],
                'resultsdataset',
                results_storage_name,
                local_cache_directory,
                suffix_filter='.h5',
            )

            for file_resource_id in results['file_resources']:
                file_resource = tantalus_api.get('file_resource', id=file_resource_id)

                if file_resource['filename'].endswith(h5_suffix):
                    for table_name, table_key, read_cols, index_cols in table_info:
                        filepath = local_results_client.get_url(file_resource['filename'])
                        data = pd.read_hdf(filepath, table_key, columns=read_cols)
                        data.set_index(index_cols, inplace=True)
                        results_tables[table_name][ticket] = data

        if subsample is not None:
            cell_ids = (
                results_tables['hmmcopy_metrics'][ticket][['cell_id']]
                .drop_duplicates().sample(frac=subsample))
            for table_name in ('hmmcopy_reads', 'hmmcopy_metrics', 'align_metrics'):
                results_tables[table_name][ticket] = results_tables[table_name][ticket].merge(cell_ids)

    for table_name, ticket_tables in results_tables.iteritems():
        print 'concatenate', table_name
        results_data = pd.concat(ticket_tables.values(), sort=True)
        results_tables[table_name] = results_data

    return results_tables


