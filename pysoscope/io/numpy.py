import numpy as np
import multiprocessing as mp 

from pysoscope.io.utils import map_peak_indices, map_peak_indices_wrapper, DEF_TOL


def convert_maldi2numpy(maldi_data, target_peaks=[], tol=DEF_TOL,ncpu=1):
    # Convert X/Z in maldi data
    x = maldi_data["data"]["x0"]
    y = maldi_data["data"]["y0"]

    Z = np.zeros((int(np.max(x)), int(np.max(y)), len(target_peaks)))

    # Map pixels to the
    this_mz = maldi_data["data"]["peak_mz"]
    this_signal = maldi_data["data"]["peak_sig"]

    if ncpu == 1 :
        for k,(xi, yi) in enumerate(zip(x,y)):
            ix, map_ix = map_peak_indices(this_mz[k], target_peaks=target_peaks, tol=tol)
            Z[int(xi-1), int(yi-1), map_ix] = this_signal[k][ix]
    else:
        args = [( [mz,], {"target_peaks":target_peaks,"tol":tol}) for mz in this_mz]
        with mp.Pool(ncpu) as pool:
            mapped_indices = pool.map(map_peak_indices_wrapper, args)
            
        for k,(xi,yi) in enumerate(zip(x,y)):
            ix, map_ix = mapped_indices[k]
            Z[int(xi-1), int(yi-1), map_ix] = this_signal[k][ix]

    return Z


def save_numpy_data(fname, **kwargs):
    np.savez(fname, **kwargs)
