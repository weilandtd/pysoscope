import numpy as np 
import multiprocessing as mp

from pysoscope.io.mat import load_maldi_file_as_dict as load_maldi_file
from pysoscope.untargeted.sample import subsample_maldi_datasets_naive
from pysoscope.untargeted.pick_peaks import get_unique_peaks
from pysoscope.untargeted.background import get_background_signature
from pysoscope.io.numpy import convert_maldi2numpy, save_numpy_data


if __name__ == "__main__":
    # Parameters and files 


    tol = 15e-6 # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names = ["EC1", "EC2", "EC3",]

    path = "./../data/"

    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
    
    # Subsample peaks and get unique peakls
    subsample = subsample_maldi_datasets_naive(maldi_datasets)

    # Get unique sample peaks
    untargeted_peaks = get_unique_peaks(subsample, tol=tol, min_samples=1)

    sample_signature_peaks = [124.007, ]

    background_sig = get_background_signature(maldi_datasets, untargeted_peaks,  sample_signature_peaks,
                                              sig_quantile=0.75,background_quantile=0.25, tol=15e-6)
    
    # Convert to numpy matrix and save
    matrix_maldi = [convert_maldi2numpy(md, target_peaks=untargeted_peaks, tol=tol, ncpu=ncpu) for md in maldi_datasets]
    [save_numpy_data("/../output/"+name+"_15ppm_min5samples_301123.npz",
                     img_data=Z, peaks=untargeted_peaks, background=background_sig) for name,Z in zip(file_names,matrix_maldi)]

    print("Done")