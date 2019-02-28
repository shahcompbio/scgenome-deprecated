import collections
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from datamanagement.miscellaneous.hdf5helper import read_python2_hdf5_dataframe


standard_hmmcopy_reads_cols = [
    'chr',
    'start',
    'end',
    'cell_id',
    'gc',
    'reads',
    'copy',
    'state',
]


categorical_cols = [
    'cell_id',
    'chr',
    'sample_id',
    'library_id',
]


def import_cn_data(
        tickets,
        local_cache_directory,
        sample_ids=None,
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
                        data = read_python2_hdf5_dataframe(filepath, table_key)
                        if read_cols is not None:
                            data = data[read_cols]
                        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
                        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]
                        if sample_ids is not None:
                            data = data[data['sample_id'].isin(sample_ids)]
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
    for table_name, ticket_tables in results_tables.items():
        for ticket, table in ticket_tables.items():
            for col in categorical_cols:
                if col not in table:
                    continue
                col_categories[col].update(table[col].cat.categories.values)

    # Create a pandas index for each set of categories
    for col, categories in col_categories.items():
        col_categories[col] = pd.Index(categories)

    # Set all categorical columns as having the same set of categories
    for table_name, ticket_tables in results_tables.items():
        for ticket, table in ticket_tables.items():
            for col in categorical_cols:
                if col not in table:
                    continue
                prev_values = table[col].astype(str).values
                table[col] = table[col].cat.set_categories(col_categories[col])
                assert not table[col].isnull().any()
                assert (prev_values == table[col].astype(str).values).all()

    # Concatenate tables, checking that categorical columns are maintained
    for table_name, ticket_tables in results_tables.items():
        print('concatenate', table_name)
        results_data = pd.concat(list(ticket_tables.values()), sort=True)
        for col in categorical_cols:
            if col not in results_data:
                continue
            assert isinstance(results_data[col].dtype, pd.api.types.CategoricalDtype)
        results_tables[table_name] = results_data

    return results_tables


