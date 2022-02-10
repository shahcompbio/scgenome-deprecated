import os

import pandas as pd
import scgenome.loaders.utils
import yaml
import scgenome.loaders.align
import scgenome.loaders.annotation
import scgenome.loaders.hmmcopy


def load_qc_files(
        hmmcopy_reads,
        hmmcopy_segs,
        hmmcopy_metrics,
        alignment_metrics,
        gc_metrics,
        annotation_metrics,
        sample_ids=None,
        additional_hmmcopy_reads_cols=None,
    ):

    results_tables = scgenome.loaders.align.load_alignment_files(
        alignment_metrics,
        gc_metrics=gc_metrics
    )

    hmmcopy_results_tables = scgenome.loaders.hmmcopy.load_hmmcopy_files(
        hmmcopy_reads,
        hmmcopy_segs, hmmcopy_metrics,
        additional_reads_cols=additional_hmmcopy_reads_cols
    )
    results_tables.update(hmmcopy_results_tables)

    annotation_results_tables = scgenome.loaders.annotation.load_annotation_files(annotation_metrics)
    results_tables.update(annotation_results_tables)

    if sample_ids is not None:
        results_tables = _sample_id_filter(results_tables, sample_ids)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_qc_results(
        alignment_results_dir,
        hmmcopy_results_dir,
        annotation_results_dir=None,
        sample_ids=None,
        additional_hmmcopy_reads_cols=None,
):
    """ Load qc data (align, hmmcopy)
    
    Args:
        alignment_results_dir (str): alignment results directory to load from.
        hmmcopy_results_dir (str): hmmcopy results directory to load from.

    KwArgs:
        annotation_results_dir (str): annotation results directory to load from. Defaults to None.
        sample_ids (list of str, optional): Set of sample ids to filter for. Defaults to None.
        additional_hmmcopy_reads_cols (list of str, optional): Additional columns to obtain from the reads table. Defaults to None.
    """

    alignment_results_dir = scgenome.loaders.utils.find_results_directory(
        alignment_results_dir, 'alignment')

    results_tables = scgenome.loaders.align.load_alignment_results(alignment_results_dir)

    hmmcopy_results_dir = scgenome.loaders.utils.find_results_directory(
        hmmcopy_results_dir, 'hmmcopy')

    hmmcopy_results_tables = scgenome.loaders.hmmcopy.load_hmmcopy_results(
        hmmcopy_results_dir,
        additional_reads_cols=additional_hmmcopy_reads_cols)

    results_tables.update(hmmcopy_results_tables)

    if annotation_results_dir is not None:
        annotation_results_dir = scgenome.loaders.utils.find_results_directory(
            annotation_results_dir, 'annotation')

        annotation_results_tables = scgenome.loaders.annotation.load_annotation_results(annotation_results_dir)

        results_tables.update(annotation_results_tables)

    if sample_ids is not None:
        results_tables = _sample_id_filter(results_tables, sample_ids)

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def _sample_id_filter(results_tables, sample_ids):
    for table_name, table_data in results_tables.items():
        results_tables[table_name] = table_data[table_data['sample_id'].isin(sample_ids)]

    return results_tables
