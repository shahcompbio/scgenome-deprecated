import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import seaborn
import numpy as np
import pandas as pd
import pylab
import sklearn.preprocessing
import scipy

import scgenome.loaders.qc
import scgenome.loaders.utils
import scgenome.loaders.hmmcopy
import scgenome.loaders.annotation
import scgenome.ds

def read_dlp_hmmcopy(hmmcopy_results_dir, annotation_results_dir, additional_hmmcopy_reads_cols=None):
    """ Read hmmcopy results from the DLP pipeline.

    Parameters
    ------
    hmmcopy_results_dir (str):
        dlp pipeline hmmcopy results directory
    annotation_results_dir (str):
        dlp pipeline annotation results directory
        
    Returns
    ------
    DLPAnnData
        An instantiated DLPAnnData Object.
    """
    
    # Load hmmcopy results
    hmmcopy_results_dir = scgenome.loaders.utils.find_results_directory(
        hmmcopy_results_dir, 'hmmcopy')

    hmmcopy_results_tables = scgenome.loaders.hmmcopy.load_hmmcopy_results(
        hmmcopy_results_dir,
        additional_reads_cols=additional_hmmcopy_reads_cols)

    results_tables = hmmcopy_results_tables

    # Load annotation results
    annotation_results_dir = scgenome.loaders.utils.find_results_directory(
        annotation_results_dir, 'annotation')

    annotation_results_tables = scgenome.loaders.annotation.load_annotation_results(annotation_results_dir)

    results_tables.update(annotation_results_tables)
    
    # Merge the metrics
    metrics = results_tables['annotation_metrics'] \
                .merge(
                    results_tables['hmmcopy_metrics'],
                    how="outer",
                    left_on='cell_id',
                    right_on='cell_id',
                    suffixes=("_leftdf", "_rightdf")
                ) 

    metrics = metrics.loc[:,~df.columns.str.endswith("_rightdf")]
    metrics = metrics.rename(
        columns = lambda x:  x.strip("_leftdf") if x.endswith("_leftdf") else x
    )
    metrics = metrics.loc[:,~metrics.columns.duplicated()]
    
    dad = DLPAnnData()
    dad.load_dlp(
        results_tables['hmmcopy_reads'],
        metrics
    )
    
    return dad
