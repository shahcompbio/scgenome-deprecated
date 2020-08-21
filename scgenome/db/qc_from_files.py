import scgenome.loaders.align
import scgenome.loaders.hmmcopy
import scgenome.loaders.qc


def get_qc_data_from_filenames(
        annotation_metrics_list, hmmcopy_reads_list, hmmcopy_segs_list,
        hmmcopy_metrics_list, alignment_metrics_list, gc_metrics_list,
        sample_ids=None, additional_hmmcopy_reads_cols=None
):
    results_tables = {}

    data = zip(annotation_metrics_list, hmmcopy_reads_list, hmmcopy_segs_list,
               hmmcopy_metrics_list, alignment_metrics_list, gc_metrics_list
               )

    for ann_metrics, hmm_reads, hmm_segs, hmm_metrics, align_metrics, gc_metrics in data:
        qc_results = scgenome.loaders.qc.load_qc_data_from_files(
            hmm_reads, hmm_segs,
            hmm_metrics, align_metrics, gc_metrics,
            annotation_metrics=ann_metrics,
            sample_id=sample_ids,
            additional_hmmcopy_reads_cols=additional_hmmcopy_reads_cols
        )

        results_tables = _aggregate_results_tables(results_tables, qc_results)

    results_tables = _concat_results_tables(results_tables)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def _concat_results_tables(results_tables):
    for table_name, table_data in results_tables.items():
        results_tables[table_name] = scgenome.utils.concat_with_categories(table_data)
    return results_tables


def _aggregate_results_tables(results_tables, ticket_results):
    for table_name, table_data in ticket_results.items():
        if table_name not in results_tables:
            results_tables[table_name] = []
        results_tables[table_name].append(table_data)

    return results_tables
