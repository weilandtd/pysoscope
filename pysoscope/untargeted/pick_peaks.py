import numpy as np 

DEF_TOL = 15e-6

def get_unique_peaks(subsample, signal_to_noise_threshold=3, tol=DEF_TOL):
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
