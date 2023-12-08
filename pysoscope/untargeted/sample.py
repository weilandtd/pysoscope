import numpy as np 
from numpy.random import permutation as randperm


def subsample_maldi_data(maldi_data, N=500):   
    n_data = len(maldi_data['data']['x'])

    # Pick N subsamples from a random index permutation
    ix = randperm(range(n_data))[1:N]
    
    sample_peaks = np.concatenate([maldi_data['data']['peak_mz'][i] for i in ix])
    sample_sigs = np.concatenate([maldi_data['data']['peak_sig'][i] for i in ix])
    
    return [(p,s) for p,s in zip(sample_peaks, sample_sigs)]


def subsample_maldi_datasets_naive(list_maldi_data, N=500):
    
    # Get subsamples 
    this_subsamples = [subsample_maldi_data(md,N=N) for md in list_maldi_data]
    
    # Annotate sample id 
    this_subsamples_id = [(p,s,i) for i, this_subsample in enumerate(this_subsamples) for p,s in this_subsample]
    
    return this_subsamples_id


# TODO Make this
def subsample_maldi_data_space(maldi_data, N=500, xlim=[0,25], ylim=[0,25], measure_from="topleft"):
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


def subsample_maldi_datasets_space(list_maldi_data, N=500, xlim=[0,25], ylim=[0,25], measure_from="topleft"):
    
    # Get subsamples 
    this_subsamples = [subsample_maldi_data_space(md,N=N,xlim=xlim,ylim=ylim,measure_from=measure_from) 
                       for md in list_maldi_data]
    
    # Annotate sample id 
    this_subsamples_id = [(p,s,i) for i, this_subsample in enumerate(this_subsamples) for p,s in this_subsample]
    
    return this_subsamples_id



