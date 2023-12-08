import numpy as np 
from numpy.random import permutation as randperm

def sample_background_pixes(maldi_data, mz_list, quantile=0.25, N=5000, tol=15e-6):
    """AI is creating summary for sample_background_pixes

    Args:
        maldi_data ([type]): [description]
        mz_list ([type]): a list of known mz values for compounds that should be in the background
        quantile (float, optional): [description]. Defaults to 0.25.
        N (int, optional): [description]. Defaults to 5000.
        tol ([type], optional): [description]. Defaults to 15e-6.
    """

    # Draw sample N    
    n_data = len(maldi_data['data']['x'])

    # Pick N subsamples from a random index permutation
    ix = randperm(range(n_data))[1:N]
    
    sample_peaks = [maldi_data['data']['peak_mz'][i] for i in ix]
    sample_sigs = [maldi_data['data']['peak_sig'][i] for i in ix]
    
    # Get mz intensities 
    mz_intensities = [[] for i in range(len(mz_list))]
    for mz_ix, mz in enumerate(mz_list):
        for sig_sample, mz_sample in zip(sample_sigs,sample_peaks):
            if any(np.abs(mz - mz_sample) <= tol*mz):
                k = np.where(np.abs(mz - mz_sample) <= tol*mz)[0]
                for i in k:
                    mz_intensities[mz_ix].append(sig_sample[i])    

    # Get quantile
    quantiles = [np.quantile(mz_intensities[i], quantile) for i in range(len(mz_list))]

    # Get sample where the respective mz intensi are everywhere below the quantile
    background_ix = [] 
    for i,(signals,peaks) in enumerate(zip(sample_sigs,sample_peaks)):
        mz_test = [False,] * len(mz_list)
        for p,s in zip(peaks,signals):
            for mz_ix, mz in enumerate(mz_list):
                if abs(mz - p) <= tol*mz and s <= quantiles[mz_ix]:
                    mz_test[mz_ix] = True

        if all(mz_test):
            background_ix.append(ix[i])

    # return background samples
    return background_ix



def get_background_signature(list_maldi_data, picked_peaks, sample_signature_peaks, 
                             sig_quantile=0.75,background_quantile=0.25, tol=15e-6 ,N=5000,):
    """AI is creating summary for get_background_signature

    Args:
        samples ([type]): List of indices of samples with background peaks 
        mz_list ([type]): A list of picked peaks 
        quantile (float, optional): [description]. Defaults to 0.75.
    """
    all_sample_sigs = []
    all_sample_peaks = []
    for maldi_data in list_maldi_data:
        background_ix = sample_background_pixes(maldi_data, sample_signature_peaks, quantile=background_quantile, N=N,tol=tol)

        sample_background_sigs = [maldi_data['data']['peak_sig'][i] for i in background_ix]
        sample_background_peaks = [maldi_data['data']['peak_mz'][i] for i in background_ix]

        all_sample_sigs.extend(sample_background_sigs)
        all_sample_peaks.extend(sample_background_peaks)
    

    # Get mz intensities 
    mz_intensities = [[] for i in range(len(picked_peaks))]

    for mz_ix, mz in enumerate(picked_peaks):
        for sig_sample,mz_sample in zip(all_sample_sigs,all_sample_peaks):
            if any(np.abs(mz - mz_sample) <= tol*mz):
                peaks_close_to_mz = np.where(np.abs(mz - mz_sample) <= tol*mz)[0]
                for peak_ix in peaks_close_to_mz:
                    mz_intensities[mz_ix].append(sig_sample[peak_ix])    
   
    # Get quantile
    quantiles = [np.quantile(mz_intensities[i], sig_quantile) for i in range(len(picked_peaks))]

    return quantiles
            


    