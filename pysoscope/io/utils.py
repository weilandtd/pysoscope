import numpy as np

DEF_TOL = 15e-6


def map_peak_indices(peaks, target_peaks=[], tol=DEF_TOL):
    """AI is creating summary for map_peak_indices

    Args:
        peaks ([list]): List of peak m/z in at a position x,y
        target_peaks (list, optional):  List of target peaks m/z
        tol (float, optional): Tolerance to match peaks. Defaults to DEF_TOL=15ppm.

    Returns:
        [type]: [description]
    """
    target_ix = []
    raw_ix = []
    for i in range(len(peaks)):
        # 
        this_ix = np.where(np.isclose(target_peaks, peaks[i], atol=tol*peaks[i]))[0]
        if len(this_ix) > 0:
            target_ix.append(this_ix[0])
            raw_ix.append(i)

    return raw_ix, target_ix


def map_peak_indices_wrapper(arg):
    args, kwargs = arg
    return map_peak_indices(*args, **kwargs)