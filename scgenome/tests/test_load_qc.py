import os
import sys
import click
import logging
import collections
import numpy as np

from scgenome.loaders.qc import load_qc_results


LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


dtypes_check = {
    'align_metrics': {
        'cell_id': 'category',
        'total_mapped_reads': 'int64',
        'library_id': 'category',
    },
    'hmmcopy_reads': {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'reads': 'int64',
        'gc': 'float64',
        'copy': 'float64',
        'state': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_segs': {
        'chr': 'category',
        'start': 'int64',
        'state': 'int64',
        'median': 'float64',
        'multiplier': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_metrics': {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'order': 'int64',
        'total_mapped_reads_hmmcopy': 'int64',
        'experimental_condition': 'object',
        'mean_copy': 'float64',
        'state_mode': 'int64',
    },
    'annotation_metrics': {
        'autocorrelation_hmmcopy': 'float64',
        'breakpoints': 'int64',
        'cell_call': 'str',
        'cell_id': 'str',
        'column': 'int64',
        'coverage_breadth': 'float64',
        'coverage_depth': 'float64',
        'cv_hmmcopy': 'float64',
        'empty_bins_hmmcopy': 'int64',
        'estimated_library_size': 'int64',
        'experimental_condition': 'str',
        'fastqscreen_grch37_multihit': 'int64',
        'fastqscreen_grch37': 'int64',
        'fastqscreen_mm10_multihit': 'int64',
        'fastqscreen_mm10': 'int64',
        'fastqscreen_nohit': 'int64',
        'fastqscreen_salmon_multihit': 'int64',
        'fastqscreen_salmon': 'int64',
        'grch37_multihit': 'int64',
        'grch37': 'int64',
        'img_col': 'int64',
        'index_i5': 'str',
        'index_i7': 'str',
        'index_sequence': 'str',
        'is_contaminated': 'bool',
        'is_s_phase_prob': 'float64',
        'is_s_phase': 'bool',
        'log_likelihood': 'float64',
        'mad_hmmcopy': 'float64',
        'mad_neutral_state': 'float64',
        'MBRSI_dispersion_non_integerness': 'float64',
        'MBRSM_dispersion': 'float64',
        'mean_copy': 'float64',
        'mean_hmmcopy_reads_per_bin': 'float64',
        'mean_insert_size': 'float64',
        'mean_state_mads': 'float64',
        'mean_state_vars': 'float64',
        'median_hmmcopy_reads_per_bin': 'float64',
        'median_insert_size': 'float64',
        'mm10_multihit': 'int64',
        'mm10': 'int64',
        'MSRSI_non_integerness': 'float64',
        'multiplier': 'int64',
        'nohit': 'int64',
        'order_corrupt_tree': 'float64',
        'order': 'int64',
        'paired_duplicate_reads': 'int64',
        'paired_mapped_reads': 'int64',
        'percent_duplicate_reads': 'float64',
        'primer_i5': 'str',
        'primer_i7': 'str',
        'quality': 'float64',
        'row': 'int64',
        'salmon_multihit': 'int64',
        'salmon': 'int64',
        'sample_id': 'str',
        'sample_type': 'str',
        'scaled_halfiness': 'float64',
        'standard_deviation_insert_size': 'float64',
        'state_mode': 'int64',
        'std_hmmcopy_reads_per_bin': 'float64',
        'total_duplicate_reads': 'int64',
        'total_halfiness': 'float64',
        'total_mapped_reads_hmmcopy': 'int64',
        'total_mapped_reads': 'int64',
        'total_properly_paired': 'int64',
        'total_reads': 'int64',
        'true_multiplier': 'float64',
        'unmapped_reads': 'int64',
        'unpaired_duplicate_reads': 'int64',
        'unpaired_mapped_reads': 'int64',
    },
    'gc_metrics': {
    },
}


def test_qc_data(results_tables):
    failed = False
    for table_name, table_data in results_tables.items():
        logging.info(f'table {table_name} has size {len(table_data)}')

        for column_name, dtype_name in dtypes_check[table_name].items():
            if column_name not in results_tables[table_name]:
                logging.warning(f'{column_name} not in table {table_name}')
                continue

            column_dtype = str(results_tables[table_name][column_name].dtype)

            expected_dtype_names = (dtype_name,)
            if dtype_name == 'str':
                expected_dtype_names = ('object', 'category')
            elif dtype_name == 'int64':
                expected_dtype_names = ('int64', 'Int64')

            if not column_dtype in expected_dtype_names:
                logging.error(f'{column_name} in {table_name} has dtype {column_dtype} not any of {expected_dtype_names}')
                failed = True

            if dtype_name == 'str' and not results_tables[table_name][column_name].apply(lambda a: isinstance(a, str)).all():
                logging.error(f'{column_name} in {table_name} is expected to be str')
                failed = True

    if failed:
        raise Exception('one or more tests failed')


def test_load_local_qc_data(results_dir):
    results_tables = load_qc_data(results_dir)
    test_qc_data(results_tables)
    logging.info(f'successfully loaded results from {results_dir}')


def test_load_stored_qc_data(tantalus_api, ticket_id, local_cache_directory=None, local_storage_name=None):
    if local_cache_directory is not None and local_storage_name is not None:
        raise ValueError('local_cache_directory and local_storage_name are mutually exclusive')

    if local_cache_directory is None and local_storage_name is None:
        raise ValueError('require one of local_cache_directory and local_storage_name')

    if local_cache_directory is not None:
        cache_qc_results(ticket_id, local_cache_directory)
        local_results_directory = local_cache_directory

    elif local_storage_name is not None:
        local_results_directory = tantalus_api.get('storage', name=local_storage_name)['storage_directory']

    ticket_directory = os.path.join(local_results_directory, ticket_id)

    test_load_local_qc_data(ticket_directory)


@click.group()
def cli():
    pass


@cli.command()
@click.argument('ticket_id')
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
def test_cached_single_ticket(ticket_id, local_cache_directory=None, local_storage_name=None):
    tantalus_api = dbclients.tantalus.TantalusApi()

    test_load_stored_qc_data(
        tantalus_api,
        ticket_id,
        local_cache_directory=local_cache_directory,
        local_storage_name=local_storage_name,
    )


@cli.command()
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
@click.option('--one_of_each')
def test_cached_multi_ticket(local_cache_directory=None, local_storage_name=None, one_of_each=False):
    tantalus_api = dbclients.tantalus.TantalusApi()

    hmmcopy_analyses = tantalus_api.list('analysis', analysis_type__name='hmmcopy')

    version_tickets = collections.defaultdict(list)
    for analysis in hmmcopy_analyses:
        version_tickets[analysis['version']].append(analysis['jira_ticket'])
    
    ticket_ids = []
    for version in version_tickets:
        if one_of_each:
            ticket_ids.append(version_tickets[version][-1])
        else:
            ticket_ids.extend(version_tickets[version])

    for ticket_id in ticket_ids:
        logging.info(ticket_id)

        ticket_results = list(tantalus_api.list('resultsdataset', analysis__jira_ticket=ticket_id))
        if len(ticket_results) == 0:
            logging.error(f'ticket {ticket_id} has no associated results')
            continue

        try:
            test_load_stored_qc_data(
                tantalus_api,
                ticket_id,
                local_cache_directory=local_cache_directory,
                local_storage_name=local_storage_name,
            )
        except KeyboardInterrupt:
            raise
        except:
            logging.exception(f"{ticket_id} failed")
        else:
            logging.exception(f"{ticket_id} succeeded")


@cli.command()
@click.argument('results_directory')
def test_local_results(results_directory):
    test_load_local_qc_data(results_directory)


@cli.command()
@click.argument('alignment_results_dir')
@click.argument('hmmcopy_results_dir')
@click.option('--annotation_results_dir')
def test_results(alignment_results_dir, hmmcopy_results_dir, annotation_results_dir=None):
    results_tables = load_qc_results(
        alignment_results_dir,
        hmmcopy_results_dir,
        annotation_results_dir=annotation_results_dir,
        sample_ids=None,
        additional_hmmcopy_reads_cols=None,
    )

    test_qc_data(results_tables)

    logging.info(f'successfully loaded results from {alignment_results_dir}, {hmmcopy_results_dir}')


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    cli()
