import anndata

import numpy as np
import pandas as pd 
import multiprocessing as mp 

#from tqdm import tqdm

from pysoscope.io.utils import map_peak_indices, map_peak_indices_wrapper, DEF_TOL


def convert_maldi_image_to_anndata(maldi_data, target_peaks=[], tol=DEF_TOL, ncpu=1, metadata=None):
    """AI is creating summary for convert_maldi2numpy

    Args:
        maldi_data ([type]): dictionary format a maldi image
        target_peaks (list, optional): A list of m/z of peaks to pick. Defaults to [].
        tol ([type], optional): Tolerance to differentiate peaks. Defaults to DEF_TOL.
        ncpu (int, optional): number of cpus. Defaults to 1.
        metadata ([dict], optional): Dictionary containing the annotation of the experimental conditions. Defaults to None.

    Returns:
        [type]: [description]
    """

    # Convert X/Z in maldi data
    x = maldi_data["data"]["x0"]
    y = maldi_data["data"]["y0"]

    data = np.zeros( (len(x), len(target_peaks)) )

    # Map pixels to the
    this_mz = maldi_data["data"]["peak_mz"]
    this_signal = maldi_data["data"]["peak_sig"]

    if ncpu == 1 :
        for k,(xi, yi) in enumerate(zip(x,y)):
            ix, map_ix = map_peak_indices(this_mz[k], target_peaks=target_peaks, tol=tol)
            data[k, map_ix] = this_signal[k][ix]
    else:
        args = [( [mz,], {"target_peaks":target_peaks,"tol":tol}) for mz in this_mz]
        with mp.Pool(ncpu) as pool:
            mapped_indices = pool.map(map_peak_indices_wrapper, args)
            
        for k,(xi,yi) in enumerate(zip(x,y)):
            ix, map_ix = mapped_indices[k]
            data[k, map_ix] = this_signal[k][ix]

    # Make sure data is numerical
    data = data.astype(np.float64)

    if metadata is None:
        obs = dict(x=x, y=y)
        var = pd.DataFrame(data = target_peaks, 
                           index=[str(round(p,4)) for p in target_peaks],
                           columns=["mz",])
        adata = anndata.AnnData(X=data, obs=obs,var=var)
    else:
        obs = dict(x=x, y=y, **metadata)
        var = pd.DataFrame(data = target_peaks, 
                           index=[str(round(p,4)) for p in target_peaks],
                           columns=["mz",])
        adata = anndata.AnnData(X=data, obs=obs, var=var)

    return adata


def join_anndata(anndata_list):
    """AI is creating summary for join_anndata

    Args:
        anndata_list ([type]): [description]

    Returns:
        [type]: [description]
    """
    return anndata.AnnData.concatenate(*anndata_list, join="outer" )
