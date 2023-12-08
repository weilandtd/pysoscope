import numpy as np 
import multiprocessing as mp

from pysoscope.io.mat import load_maldi_file_as_dict as load_maldi_file
from pysoscope.untargeted.sample import subsample_maldi_datasets_naive
from pysoscope.untargeted.pick_peaks import get_unique_peaks
from pysoscope.io.numpy import convert_maldi2numpy, save_numpy_data

if __name__ == "__main__":
    # Parameters and files 

    tol = 15e-6  # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names =["EC1_1", "EC2_1", "EC2_2", "ED1_1","ED1_2", "ED2_1","EDL2_2","ED3_1", "ED3_2", "SC1_1" ,"SC2_1","SD1_1", "SD1_2", "SD2_1", "SD2_2", "SD3_1","SD3_2"]
#["20211025-FBP-20um-normal-infusion_new","20211025-FBP-20um-pulse-chase-1"]
    path = "/scratch/gpfs/cm7897/raw mat files/"

    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
    
    # Subsample peaks and get unique peakls
    subsample = subsample_maldi_datasets_naive(maldi_datasets)

    # Get unique sample peaks
    real_unique_peaks_filtered = get_unique_peaks(subsample, tol=tol, min_samples=3)
    
    # Convert to numpy matrix and save
    matrix_maldi = [convert_maldi2numpy(md, target_peaks=real_unique_peaks_filtered, tol=tol, ncpu=ncpu) for md in maldi_datasets]
    [save_numpy_data("/scratch/gpfs/cm7897/output/"+name+"_15ppm_min_samples_5.npz",img_data=Z, peaks=real_unique_peaks_filtered) for name,Z in zip(file_names,matrix_maldi)]

    print("Done")
