import numpy as np 
from numpy.random import permutation as randperm

import h5py
import multiprocessing as mp

tol = 15e-6


def load_maldi_file(filepath,):
    maldi_data = {"data":{} ,"fname":[], "res":[]}
    f = h5py.File(filepath)
    # Load meta data
    msi = f['msi']
    data = msi['data']
    maldi_data['fname'] = np.array(msi['fname'])
    maldi_data['res'] = np.array(msi['res'])
    # Read the data
    for k, v in data.items():
        if not ( k in ["peak_mz"], ["peak_sig"]):
            maldi_data["data"][k] = [np.array(f[vi[0]])[0][0] for vi in v ]
        else:
            maldi_data["data"][k] = [np.array(f[vi[0]])[0] for vi in v ]
 
    # Convert X/Z in maldi data 
    _x = maldi_data["data"]["x"]
    _y = maldi_data["data"]["y"]

    min_x = np.min(_x)
    min_y = np.min(_y)

    x = ( _x - min_x) + 1 ;
    y = ( _y - min_y) + 1 ;

    print(x[:,0])
    maldi_data["data"]["x0"] = x[:,0]
    maldi_data["data"]["y0"] = y[:,0]
    
    return maldi_data

def subsample_maldi_file(maldi_data, N=500):   
    n_data = len(maldi_data['data']['x'])

    # Pick N subsamples from a random index permutation
    ix = randperm(range(n_data))[1:N]
    
    sample_peaks = np.concatenate([maldi_data['data']['peak_mz'][i] for i in ix])
    sample_sigs = np.concatenate([maldi_data['data']['peak_sig'][i] for i in ix])
    
    return [(p,s) for p,s in zip(sample_peaks, sample_sigs)]

def subsample_maldi_file_space(maldi_data, N=500, xlim=[0,25], ylim=[0,25], measure_from="topleft"):
    x = maldi_data["data"]["x0"]
    y = maldi_data["data"]["y0"]
    
    max_x = np.max(x)
    max_y = np.max(y)

    if measure_from == "topleft":
        xlim = [max_x - xlim[1], max_x - xlim[0]]
        ylim = [max_y - ylim[1], max_y - ylim[0]]
    elif measure_from =="center":
        xlim = [max_x/2 - xlim[1], max_x/2 - xlim[0]]
        ylim = [max_y/2 - ylim[1], max_y/2 - ylim[0]]
    else:
        raise NotImplemented("KWD arg not implemented!")

    #print(xlim,ylim, max_x, max_y)
    region_ix = [] 
    for i, (xi,yi) in enumerate(zip(x,y)):
        if (xi >= xlim[0]) and (xi <= xlim[1]) and (yi >= ylim[0]) and (yi <= ylim[1]):
            region_ix.append(i)

    #print(len(region_ix))
    # Pick subsample 
    if len(region_ix) < N:
        raise ValueError("More samples asked than pixels in defined region!")
        
    # Pick N subsamples from a random index permutation
    ix = randperm(region_ix)[1:N]

    sample_peaks = np.concatenate([maldi_data['data']['peak_mz'][i] for i in ix])
    sample_sigs = np.concatenate([maldi_data['data']['peak_sig'][i] for i in ix])
    
    return [(p,s) for p,s in zip(sample_peaks, sample_sigs)]


def get_unique_peaks(subsample, signal_to_noise_threshold=3, tol = 15e-6):
    # Sort peaks subsamples by peaks
    peaks, signals = zip( *sorted(subsample, key=lambda pair: pair[0],))

    unique_peaks = [peaks[0], ]
    unique_peaks_signals = [ [signals[0],], ]
    j = 0 
    for i in range(1,len(peaks)):
        if abs(peaks[i-1] - peaks[i] ) >= tol * peaks[i-1]:
            unique_peaks[-1] = np.median(peaks[j:i])
            unique_peaks.append(peaks[i])
            unique_peaks_signals.append([signals[i],] )
            j = i
        else:
            unique_peaks_signals[-1].append(signals[i])

    unique_peaks_expanded = np.concatenate([np.array([p,]*len(s)) for p,s in zip(unique_peaks,unique_peaks_signals) ] )
    unique_signals_expanded = np.concatenate([np.array(s) for s in unique_peaks_signals])
    
    cutoff = np.quantile(unique_signals_expanded, 0.8) * signal_to_noise_threshold

    unique_peaks_filtered = np.array( list(set( unique_peaks_expanded[ unique_signals_expanded > cutoff ] )) ) 

    return unique_peaks_filtered


def exclude_background_peaks(unique_peaks_sample,unique_peaks_background, tol=15e-6):
    real_unique_peaks_filtered = []
    for p in unique_peaks_sample:
        p_in_q = [] 
        for q in unique_peaks_background:
            if abs(p-q) > tol* p:
                p_in_q.append(False)
            else:
                p_in_q.append(True)
        if not any(p_in_q):
            real_unique_peaks_filtered.append(p)
            
    return real_unique_peaks_filtered


def map_peak_indices(peaks, unique_peaks=[], tol=tol):
    mapped_ix = []
    ix = []
    for i in range(len(peaks)):
        this_ix = np.where(np.isclose(unique_peaks, peaks[i], atol=tol*peaks[i]))[0]
        if len(this_ix) > 0:
            mapped_ix.append(this_ix[0])
            ix.append(i)
    return ix, mapped_ix


def map_peak_indices_wrapper(arg):
    args, kwargs = arg
    return map_peak_indices(*args, **kwargs)


def convert_maldi2matrix(maldi_data, target_peaks=[], tol=tol,ncpu=1):
    # Convert X/Z in maldi data
    x = maldi_data["data"]["x0"]
    y = maldi_data["data"]["y0"]

    Z = np.zeros((int(np.max(x)), int(np.max(y)), len(target_peaks)))

    # Map pixels to the
    this_mz = maldi_data["data"]["peak_mz"]
    this_signal = maldi_data["data"]["peak_sig"]

    if ncpu == 1 :
        for k,(xi, yi) in enumerate(zip(x,y)):
            ix, map_ix = map_peak_indices(this_mz[k], unique_peaks=target_peaks, tol=tol)
            Z[int(xi-1), int(yi-1), map_ix] = this_signal[k][ix]
    else:
        args = [( [mz,], {"unique_peaks":target_peaks,"tol":tol}) for mz in this_mz]
        with mp.Pool(ncpu) as pool:
            mapped_indices = pool.map(map_peak_indices_wrapper, args)
            
        for k,(xi,yi) in enumerate(zip(x,y)):
            ix, map_ix = mapped_indices[k]
            Z[int(xi-1), int(yi-1), map_ix] = this_signal[k][ix]

    return Z


if __name__ == "__main__":

    ncpu = 36
    file_names = ["R1_9AA", "R3_9AA", "R5_9AA", "R7_9AA",
                  "R2_1_9AA","R6_1_9AA","R7_1_9AA", "R2_9AA", "R8_9AA" ,"R4_9AA",]
    path = "./cm_data/"

    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(10) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
        

    subsample_bgd = np.concatenate([subsample_maldi_file_space(md,N=500) for md in maldi_datasets[:4]])
    subsample_sig = np.concatenate([subsample_maldi_file_space(md,ylim=[-15,15], xlim=[-15,15], N=500,measure_from="center") 
                                    for md in maldi_datasets[:4]])

    unique_peaks_background = get_unique_peaks(subsample_bgd)
    unique_peaks_sample = get_unique_peaks(subsample_sig)

    real_unique_peaks_filtered = exclude_background_peaks(unique_peaks_sample, unique_peaks_background)
    print(real_unique_peaks_filtered)
    
    matrix_maldi = [convert_maldi2matrix(md, target_peaks=real_unique_peaks_filtered, tol=tol, ncpu=ncpu) for md in maldi_datasets]

    [np.savez_compressed("./cm_data/"+name+".npz",Z) for name,Z in zip(file_names,matrix_maldi)]
    np.savez_compressed("./cm_data/peaks.npz", real_unique_peaks_filtered)

    print("Done")
    
    
    