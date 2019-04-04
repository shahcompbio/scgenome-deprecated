import collections
import logging
import pandas as pd

import dbclients.tantalus
import datamanagement.transfer_files
from dbclients.basicclient import NotFoundError


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


class HMMCopyData:
    def __init__(
            self,
            ticket_id,
            local_cache_directory,
            ploidy_solution='0',
            results_storage_name='singlecellblob_results',
    ):
        self.local_storage_directory = local_cache_directory
        self.results_storage_name = results_storage_name

        self.results_info = [
            ('hmmcopy_reads', '{}_hmmcopy', f'_multiplier{ploidy_solution}_reads.csv.gz'),
            ('hmmcopy_metrics', '{}_hmmcopy', f'_multiplier{ploidy_solution}_metrics.csv.gz'),
            ('align_metrics', '{}_align', '_alignment_metrics.csv.gz'),
        ]

        self.results_filepaths = {}

        tantalus_api = dbclients.tantalus.TantalusApi()

        for table_name, analysis_template, csv_suffix in self.results_info:
            results = tantalus_api.get('results', name=analysis_template.format(ticket_id))

            assert table_name not in self.results_filepaths

            filepaths = datamanagement.transfer_files.cache_dataset(
                tantalus_api,
                results['id'],
                'resultsdataset',
                results_storage_name,
                local_cache_directory,
                suffix_filter=csv_suffix,
            )

            if len(filepaths) != 1:
                raise Exception(f'found {len(filepaths)} filepaths for {csv_suffix}, results {results["id"]}')

            self.results_filepaths[table_name] = filepaths[0]

    def load_cn_data(self, sample_ids=None, additional_reads_cols=None, subsample=None):
        """ Load copy number and align tables
        
        Args:
            sample_ids (list of str, optional): Set of sample ids to filter for. Defaults to None.
            additional_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
            subsample (float, optional): Proportion of the cells to downsample to. Defaults to None.
        
        Returns:
            dict: pandas.DataFrame tables keyed by table name`
        """
        results_tables = {}

        hmmcopy_reads_cols = standard_hmmcopy_reads_cols
        if additional_reads_cols is not None:
            hmmcopy_reads_cols.extend(additional_reads_cols)

        for table_name, filepath in self.results_filepaths.items():
            usecols = None
            if table_name == 'hmmcopy_reads':
                usecols = hmmcopy_reads_cols

            data = pd.read_csv(filepath, usecols=usecols)

            data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
            data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

            if sample_ids is not None:
                data = data[data['sample_id'].isin(sample_ids)]

            for col in categorical_cols:
                if col in data:
                    data[col] = pd.Categorical(data[col])

            results_tables[table_name] = data

        if subsample is not None:
            cell_ids = (
                results_tables['hmmcopy_metrics'][['cell_id']]
                .drop_duplicates().sample(frac=subsample))
            for table_name in ('hmmcopy_reads', 'hmmcopy_metrics', 'align_metrics'):
                results_tables[table_name] = results_tables[table_name].merge(cell_ids)

        # For columns that need to be categorical, create a set of categories
        col_categories = collections.defaultdict(set)
        for table_name, table in results_tables.items():
            for col in categorical_cols:
                if col not in table:
                    continue
                col_categories[col].update(table[col].cat.categories.values)

        # Create a pandas index for each set of categories
        for col, categories in col_categories.items():
            col_categories[col] = pd.Index(categories)

        # Set all categorical columns as having the same set of categories
        for table_name, table in results_tables.items():
            for col in categorical_cols:
                if col not in table:
                    continue
                prev_values = table[col].astype(str).values
                table[col] = table[col].cat.set_categories(col_categories[col])
                assert not table[col].isnull().any()
                assert (prev_values == table[col].astype(str).values).all()

        # # Concatenate tables, checking that categorical columns are maintained
        # for table_name, ticket_tables in results_tables.items():
        #     print('concatenate', table_name)
        #     results_data = pd.concat(list(ticket_tables.values()), sort=True)
        #     for col in categorical_cols:
        #         if col not in results_data:
        #             continue
        #         assert isinstance(results_data[col].dtype, pd.api.types.CategoricalDtype)
        #     results_tables[table_name] = results_data

        if 'total_mapped_reads_hmmcopy' not in results_tables['hmmcopy_metrics']:
            results_tables['hmmcopy_metrics']['total_mapped_reads_hmmcopy'] = results_tables['hmmcopy_metrics']['total_mapped_reads']

        elif results_tables['hmmcopy_metrics']['total_mapped_reads_hmmcopy'].isnull().any():
            fix_read_count = results_tables['hmmcopy_metrics']['total_mapped_reads_hmmcopy'].isnull()
            results_tables['hmmcopy_metrics'].loc[fix_read_count, 'total_mapped_reads_hmmcopy'] = (
                results_tables['hmmcopy_metrics'].loc[fix_read_count, 'total_mapped_reads'])

        return results_tables

