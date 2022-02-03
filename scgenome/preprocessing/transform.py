import numpy as np

from numpy import ndarray


def fill_missing(data: ndarray) -> ndarray:
    """ Fill missing values with means of non-nan values across rows.

    Parameters
    ----------
    data : ndarray
        data to fill

    Returns
    -------
    ndarray
        copy of data with missing entries filled 
    """    

    # Deal with missing values by assigning the mean value
    # of each bin to missing values of that bin
    data = data.copy()
    bin_means = np.nanmean(data, axis=0)
    bin_means = np.nan_to_num(bin_means, nan=0)
    bin_means = np.tile(bin_means, (data.shape[0], 1))
    data[np.where(np.isnan(data))] = bin_means[np.where(np.isnan(data))]

    return data

