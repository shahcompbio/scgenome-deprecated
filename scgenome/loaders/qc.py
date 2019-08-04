import click
import scgenome.loaders.utils
from scgenome.loaders.align import load_align_data
from scgenome.loaders.hmmcopy import load_hmmcopy_data
from scgenome.loaders.annotation import load_annotation_data


def load_cached_qc_data(
        ticket_id,
        local_cache_directory,
        sample_ids=None,
        subsample=None,
        additional_hmmcopy_reads_cols=None,
    ):
    """ Load qc data (align, hmmcopy, annotation)
    
    Args:
        ticket_id (str): jira ticket for the analyis producing the results.
        local_cache_directory (str): path to hmmcopy results.
        sample_ids (list of str, optional): Set of sample ids to filter for. Defaults to None.
        subsample (float, optional): Proportion of the cells to downsample to. Defaults to None.
        additional_hmmcopy_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    """

    ticket_results_dirs = scgenome.loaders.utils.find_results_directories(
        ticket_id, local_cache_directory)

    results_tables = load_align_data(ticket_results_dirs['align'])

    hmmcopy_results_tables = load_hmmcopy_data(
        ticket_results_dirs['hmmcopy'],
        additional_reads_cols=additional_hmmcopy_reads_cols)

    results_tables.update(hmmcopy_results_tables)

    # Annotation table overrides hmmcopy_metrics
    if 'annotation' in ticket_results_dirs:
        annotation_results_tables = load_annotation_data(ticket_results_dirs['annotation'])
        results_tables['hmmcopy_metrics'] = annotation_results_tables['annotation_metrics']

    # Optionally select specific samples
    if sample_ids is not None:
        for table_name, table_data in results_tables.items():
            table_data = table_data[table_data['sample_id'].isin(sample_ids)]

    # Optionally subsample cells
    if subsample is not None:
        cell_ids = (
            results_tables['hmmcopy_metrics'][['cell_id']]
            .drop_duplicates().sample(frac=subsample))
        for table_name in results_tables.keys():
            results_tables[table_name] = results_tables[table_name].merge(cell_ids)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


@click.command()
@click.argument('ticket_id')
@click.argument('local_cache_directory')
def test_load_cached_qc_data(ticket_id, local_cache_directory):
    load_cached_qc_data(
        ticket_id,
        local_cache_directory,
    )

    for table_name, table_data in results_tables.items():
        print(table_name)
        print(table_data.shape)
        print(table_data.head())


if __name__ == '__main__':
    test_load_cached_qc_data()
