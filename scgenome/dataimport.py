import collections
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files


standard_hmmcopy_reads_cols = [
    'chr',
    'start',
    'end',
    'cell_id',
    'gc',
    'reads',
    'copy',
    'state',
    'integer_copy_number',
]


categorical_cols = [
    'cell_id',
    'chr',
]


def import_cn_data(
        tickets,
        local_cache_directory,
        ploidy_solution='0',
        additional_reads_cols=None,
        results_storage_name='singlecellblob_results',
        subsample=None,
):
    tantalus_api = dbclients.tantalus.TantalusApi()

    local_results_client = tantalus_api.get_cache_client(local_cache_directory)

    results_tables = collections.defaultdict(dict)

    hmmcopy_reads_cols = standard_hmmcopy_reads_cols
    if additional_reads_cols is not None:
        hmmcopy_reads_cols.extend(additional_reads_cols)

    results_info = [
        ('{}_hmmcopy', '_hmmcopy.h5', [
            ('hmmcopy_reads', '/hmmcopy/reads/'+ploidy_solution, hmmcopy_reads_cols),
            ('hmmcopy_metrics', '/hmmcopy/metrics/'+ploidy_solution, None),
        ]),
        ('{}_align', '_alignment_metrics.h5', [
            ('align_metrics', '/alignment/metrics', None),
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
                    for table_name, table_key, read_cols in table_info:
                        filepath = local_results_client.get_url(file_resource['filename'])
                        data = pd.read_hdf(filepath, table_key, columns=read_cols)
                        for col in categorical_cols:
                            if col in data:
                                data[col] = pd.Categorical(data[col])
                        results_tables[table_name][ticket] = data

        if subsample is not None:
            cell_ids = (
                results_tables['hmmcopy_metrics'][ticket][['cell_id']]
                .drop_duplicates().sample(frac=subsample))
            for table_name in ('hmmcopy_reads', 'hmmcopy_metrics', 'align_metrics'):
                results_tables[table_name][ticket] = results_tables[table_name][ticket].merge(cell_ids)

    # For columns that need to be categorical, create a set of categories
    col_categories = collections.defaultdict(set)
    for table_name, ticket_tables in results_tables.iteritems():
        for ticket, table in ticket_tables.iteritems():
            for col in categorical_cols:
                if col not in table:
                    continue
                col_categories[col].update(table[col].cat.categories.values)

    # Create a pandas index for each set of categories
    for col, categories in col_categories.iteritems():
        col_categories[col] = pd.Index(categories)

    # Set all categorical columns as having the same set of categories
    for table_name, ticket_tables in results_tables.iteritems():
        for ticket, table in ticket_tables.iteritems():
            for col in categorical_cols:
                if col not in table:
                    continue
                prev_values = table[col].astype(str).values
                table[col] = table[col].cat.set_categories(col_categories[col])
                assert not table[col].isnull().any()
                assert (prev_values == table[col].astype(str).values).all()

    # Concatenate tables, checking that categorical columns are maintained
    for table_name, ticket_tables in results_tables.iteritems():
        print 'concatenate', table_name
        results_data = pd.concat(ticket_tables.values(), sort=True)
        for col in categorical_cols:
            if col not in results_data:
                continue
            assert isinstance(results_data[col].dtype, pd.api.types.CategoricalDtype)
        results_tables[table_name] = results_data

    return results_tables


