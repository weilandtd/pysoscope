import numpy as np 

DEF_TOL = 15e-6

def get_unique_peaks(subsample, signal_to_noise_threshold=3, tol=15e-6, min_samples=5):
    # Sort peaks subsamples by peaks
    peaks, signals, sample = zip( *sorted(subsample, key=lambda pair: pair[0],))

    peaks = np.array(peaks)
    signals = np.array(signals)
    sample = np.array(sample)
    
    # Combine all peaks in to one list and filter by signal to noise ratio 
    #unique_peaks_expanded = np.concatenate([np.array([p,]*len(s)) for p,s in zip(unique_peaks,unique_peaks_signals) ] )
    #unique_signals_expanded = np.concatenate([np.array(s) for s in unique_peaks_signals])
    
    cutoff = np.quantile(signals, 0.8) * signal_to_noise_threshold
    peaks_filtered = peaks[ signals > cutoff ]
    samples_filtered = sample[ signals > cutoff ]

    #peaks_filtered_sorted = sorted(peaks_filtered)
    
    # Filter peaks based on ppm 
    unique_peaks = [peaks_filtered[0], ]
    unique_peaks_samples = [ [samples_filtered[0],], ]
    j = 0 
    for i in range(1,len(peaks_filtered)):
        if abs(peaks_filtered[i-1] - peaks_filtered[i] ) >= tol * peaks_filtered[i-1]:
            unique_peaks[-1] = np.median(peaks_filtered[j:i])
            unique_peaks.append(peaks_filtered[i])
            unique_peaks_samples.append([samples_filtered[i],] )
            
            j = i
        else:
            unique_peaks_samples[-1].append(samples_filtered[i])

    unique_peaks_sample = [p for p, s in zip(unique_peaks, unique_peaks_samples) if len(s) > min_samples ]
    #unique_peaks_filtered = np.array( list(set( unique_peaks_expanded[ unique_signals_expanded > cutoff ] )) ) 

    return unique_peaks_sample


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
