

import anndata as ad
import pkg_resources

from anndata import AnnData


def OV2295_HMMCopy_reduced() -> AnnData:
    """ DLP data from the OV2295 ovarian cell lines.

    Returns
    -------
    AnnData
        HMMCopy data, reduced size
    """

    adata_filename = pkg_resources.resource_filename('scgenome', 'datasets/data/OV2295_HMMCopy_reduced.h5ad')
    return ad.read(adata_filename)


def OV_051_Medicc2_reduced() -> AnnData:
    """ DLP data from the OV2295 ovarian cell lines.

    Returns
    -------
    AnnData
        HMMCopy data, reduced size
    """

    adata_filename = pkg_resources.resource_filename('scgenome', 'datasets/data/OV_051_Medicc2_reduced.h5ad')
    return ad.read(adata_filename)
