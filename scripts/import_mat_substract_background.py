import numpy as np 
import multiprocessing as mp

from pysoscope.io.mat import load_maldi_file_as_dict as load_maldi_file
from pysoscope.untargeted.sample import subsample_maldi_datasets_space
from pysoscope.untargeted.pick_peaks import get_unique_peaks, exclude_background_peaks
from pysoscope.io.numpy import convert_maldi2numpy, save_numpy_data


if __name__ == "__main__":

    tol = 15e-6 # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names = ["EC1", "EC2", "EC3", "EC4", "ED1","ED2","ED3", "ED4", "PC1" ,"PC2","PC3", "PC4", "PD1","PD2", "PD3", "PD4", "SC1", "SC2", "SC3", "SC4", "SD1", "SD2", "SD3", "SD4"]
    
    # Define the sample that have background 
    samples_with_background = ["EC1", "ED4", "PC2", "PD3","SC1","SD4"]

    path = "/scratch/gpfs/cm7897/raw mat files/NEDC final/"


    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
        
    ## Subsample sample and backaground 
    index_samples_with_background = [file_names.index(name) for name in samples_with_background]
    subsample_bgd = subsample_maldi_datasets_space([maldi_datasets[i] for i in index_samples_with_background], N=500)
    subsample_sig = subsample_maldi_datasets_space([maldi_datasets[i] for i in index_samples_with_background], 
                                                   ylim=[-15,15], xlim=[-15,15], N=500,measure_from="center") 

    # Filter peaks 
    unique_peaks_background = get_unique_peaks(subsample_bgd, tol=tol, signal_to_noise_threshold=9) #Calculated based on mean signal in sample and bgd
    unique_peaks_sample = get_unique_peaks(subsample_sig, tol=tol, signal_to_noise_threshold=3)

    # Get unique sample peaks
    real_unique_peaks_filtered = exclude_background_peaks(unique_peaks_sample, unique_peaks_background,tol=tol)
    print(real_unique_peaks_filtered)

    # Convert to numpy matrix and save
    matrix_maldi = [convert_maldi2numpy(md, target_peaks=real_unique_peaks_filtered, tol=tol, ncpu=ncpu) for md in maldi_datasets]
    [save_numpy_data("/scratch/gpfs/cm7897/output/"+name+"_15ppm_min5samples_301123.npz",img_data=Z, peaks=real_unique_peaks_filtered) for name,Z in zip(file_names,matrix_maldi)]

    print("Done")
    
    
    
