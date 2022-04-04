import anndata as ad
import pandas as pd

from typing import List
from anndata import AnnData


def ad_concat_cells(adatas: List[AnnData]) -> AnnData:
    """ Concatenate a list of anndata by obs (cells)

    Parameters
    ----------
    adatas : List[AnnData]
        list of anndata to concatenate

    Returns
    -------
    AnnData
        concatenated anndata with all obs from the input anndatas
    """    

    var_temp = pd.concat((a.var for a in adatas)).drop_duplicates()
    adata = ad.concat(adatas, join='outer')
    adata.var = var_temp.reindex(adata.var.index)

    return adata
