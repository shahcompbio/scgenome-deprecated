import logging
import pandas as pd

import datamanagement.transfer_files
import dbclients.tantalus

import scgenome.snpdata
import scgenome.snvdata


class PseudobulkData:
    def __init__(
            self,
            ticket_id,
            local_cache_directory,
            results_storage_name='singlecellblob_results',
    ):
        tantalus_api = dbclients.tantalus.TantalusApi()

        self.ticket_id = ticket_id
        self.results = search_pseudobulk_results(tantalus_api, ticket_id)
        self.sample_libraries = get_pseudobulk_sample_libraries(tantalus_api, ticket_id)

        self.dataset_filepaths = datamanagement.transfer_files.cache_dataset(
            tantalus_api,
            self.results['id'],
            'resultsdataset',
            results_storage_name,
            local_cache_directory,
        )

    def get_pseudobulk_files(self, suffix):
        """ Iterate through sample / library specific files with a given suffix.
        
        Args:
            dataset_filepaths (list of str): absolute paths to pseudobulk results files
            sample_libraries (list of tuples): sample id, library id tuples
            suffix (str): suffix of filepaths to return
        
        Raises:
            ValueError: found 0 or multiple files
        """
        for sample_id, library_id in self.sample_libraries:
            sample_lib_suffix = f'{sample_id}_{library_id}_{suffix}'
            sample_lib_filepaths = list(filter(lambda a: a.endswith(sample_lib_suffix), self.dataset_filepaths))

            if len(sample_lib_filepaths) != 1:
                raise ValueError(f'found {len(sample_lib_filepaths)} {suffix} files for {sample_id}, {library_id}, {self.ticket_id}')

            sample_lib_filepath = sample_lib_filepaths[0]

            yield sample_id, library_id, sample_lib_filepath

    def load_snv_count_data(self, positions):
        """ Load per cell SNV count data
        
        Args:
            positions (pandas.DataFrame): restrict to the specified positions
        
        Returns:
            pandas.DataFrame: SNV alt and ref counts per cell
        """
        snv_count_data = []

        for sample_id, library_id, filepath in self.get_pseudobulk_files('snv_union_counts.csv.gz'):
            logging.info('Loading snv counts from {}'.format(filepath))

            data = pd.read_csv(
                filepath,
                dtype={
                    'chrom': 'category',
                    'ref': 'category',
                    'alt': 'category',
                    'cell_id': 'category',
                })

            logging.info('Loaded snv counts table with shape {}'.format(data.shape))

            scgenome.utils.union_categories(
                data, positions,
                cols=['chrom', 'ref', 'alt'])
            data = data.merge(positions, how='inner')

            logging.info('Filtered snv counts table to shape {}'.format(data.shape))

            snv_count_data.append(data)

        snv_count_data = scgenome.utils.concat_with_categories(snv_count_data, ignore_index=True)

        logging.info('Loaded all snv counts tables with shape {}'.format(snv_count_data.shape))

        return snv_count_data

    def load_snv_annotation_data(self, table_name):
        """ Load SNV annotation data
        
        Args:
            table_name (str): name of annotation table to load.
        
        Returns:
            pandas.DataFrame: SNVs annotation data per sample / library
        """
        
        snv_data = []
        for sample_id, library_id, filepath in self.get_pseudobulk_files(f'snv_{table_name}.csv.gz'):
            logging.info(f'Loading snv {table_name} annotations from {filepath}')

            data = pd.read_csv(
                filepath,
                dtype={
                    'chrom': 'category',
                    'ref': 'category',
                    'alt': 'category',
                    'cell_id': 'category',
                    'effect': 'category',
                    'effect_impact': 'category',
                    'functional_class': 'category',
                    'codon_change': 'category',
                    'amino_acid_change': 'category',
                    'gene_name': 'category',
                    'transcript_biotype': 'category',
                    'gene_coding': 'category',
                    'transcript_id': 'category',
                    'genotype': 'category',
                })

            snv_data.append(data)

        snv_data = scgenome.utils.concat_with_categories(snv_data, ignore_index=True)

        return snv_data

    def load_haplotype_allele_counts(self):
        """ Load the haplotype allele count data from the pseudobulk results paths
        
        Returns:
            pandas.DataFrame: Haplotype allele counts per cell
        """
        allele_counts = []

        for sample_id, library_id, filepath in self.get_pseudobulk_files('allele_counts.csv'):
            logging.info('Loading haplotype allele counts from {}'.format(filepath))

            # HACK: temp fix until the pipeline outputs csv headers
            names = None
            firstline = open(filepath).readline()
            if 'chromosome' not in firstline:
                names = [
                    'start',
                    'end',
                    'hap_label',
                    'allele_id',
                    'readcount',
                    'chromosome',
                    'cell_id',
                ]

            data = pd.read_csv(
                filepath,
                names=names,
                dtype={
                    'chromosome': 'category',
                    'cell_id': 'category',
                })

            logging.info('Loaded haplotype allele counts table with shape {}'.format(data.shape))

            allele_counts.append(data)

        allele_counts = scgenome.utils.concat_with_categories(allele_counts, ignore_index=True)

        logging.info('Loaded all haplotype allele counts table with shape {}'.format(allele_counts.shape))

        return allele_counts

    def load_allele_data(
            self,
            clusters,
            cn_bin_size=500000,
    ):
        """ Load cluster specific haplotype allele counts
        
        Args:
            clusters (pandas.DataFrame): cell clusters, columns cell_id, cluster_id
            cn_bin_size (int, optional): HMMCopy bin size. Defaults to 500000.
        
        Returns:
            pandas.DataFrame: haplotype allele specific counts per cluster
        """
        allele_data = self.load_haplotype_allele_counts()

        allele_data = scgenome.snpdata.calculate_cluster_allele_counts(
            allele_data, clusters, cn_bin_size)

        return allele_data

    def load_breakpoint_data(self):
        """ Load breakpoint data from a pseudobulk run.
        
        Args:
            dataset_filepaths (list of str): absolute paths to pseudobulk results files
            sample_libraries (list of tuples): sample id, library id tuples
        """
        
        breakpoint_data = []

        for sample_id, library_id, filepath in self.get_pseudobulk_files('destruct.csv.gz'):
            data = pd.read_csv(filepath, sep='\t')
            data['library_id'] = library_id
            data['sample_id'] = sample_id
            breakpoint_data.append(data)

        if len(breakpoint_data) == 0:
            return pd.DataFrame(), pd.DataFrame()

        breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

        breakpoint_count_data = []

        for sample_id, library_id, filepath in self.get_pseudobulk_files('cell_counts_destruct.csv.gz'):
            data = pd.read_csv(filepath, names=['cluster_id', 'cell_id', 'read_count'])
            data['library_id'] = library_id
            data['sample_id'] = sample_id
            breakpoint_count_data.append(data)

        breakpoint_count_data = pd.concat(breakpoint_count_data, ignore_index=True)
        breakpoint_count_data = breakpoint_count_data.rename(columns={'cluster_id': 'prediction_id'})

        # KLUDGE: normal reads are not filtered properly, filter by their prefix
        breakpoint_count_data = breakpoint_count_data[~breakpoint_count_data['cell_id'].str.startswith('HS')]

        return breakpoint_data, breakpoint_count_data


def search_pseudobulk_results(
        tantalus_api,
        ticket_id,
):
    """ Search for pseudobulk results and analysis for a given ticket.
    """
    analysis = tantalus_api.get(
        'analysis',
        jira_ticket=ticket_id,
    )

    results = tantalus_api.get(
        'resultsdataset',
        analysis=analysis['id'],
    )

    return results


def get_pseudobulk_sample_libraries(
    tantalus_api,
    ticket_id,
):
    """ Get a list of library / sample pairs for a pseudobulk library.
    """
    analysis = tantalus_api.get(
        'analysis',
        jira_ticket=ticket_id,
    )

    sample_libraries = []
    for dataset_id in analysis['input_datasets']:
        dataset = tantalus_api.get('sequencedataset', id=dataset_id)

        sample_id = dataset['sample']['sample_id']
        library_id = dataset['library']['library_id']

        if library_id == analysis['args']['matched_normal_library']:
            assert sample_id == analysis['args']['matched_normal_sample']
            continue

        sample_libraries.append((sample_id, library_id))

    return sample_libraries


