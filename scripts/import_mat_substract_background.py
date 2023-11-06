import numpy as np 
import multiprocessing as mp

from pysoscope.io.mat import load_maldi_file_as_dict as load_maldi_file
from pysoscope.untargeted.sample import subsample_maldi_file_space
from pysoscope.untargeted.pick_peaks import get_unique_peaks, exclude_background_peaks
from pysoscope.io.numpy import convert_maldi2numpy, save_numpy_data


if __name__ == "__main__":

    tol = 15e-6 # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names = ["R1_9AA", "R3_9AA", "R5_9AA", "R7_9AA", "R2_1_9AA","R6_9AA","R7_1_9AA", "R2_9AA", "R8_9AA" ,"R4_9AA",]
    
    # Define the sample that have background 
    samples_with_background = ["R1_9AA", "R3_9AA", "R5_9AA", "R7_9AA",]

    path = "/scratch/gpfs/cm7897/raw mat files/"


    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
        
    ## Subsample sample and backaground 
    index_samples_with_background = [file_names.index(name) for name in samples_with_background]
    subsample_bgd = np.concatenate([subsample_maldi_file_space(maldi_datasets[i],N=500) 
                                    for i in index_samples_with_background])
    
    subsample_sig = np.concatenate([subsample_maldi_file_space(maldi_datasets[i],
			            ylim=[-15,15], xlim=[-15,15], N=500,measure_from="center") 
                                    for i in index_samples_with_background])

    # Filter peaks 
    unique_peaks_background = get_unique_peaks(subsample_bgd, tol=tol)
    unique_peaks_sample = get_unique_peaks(subsample_sig, tol=tol)

    # Get unique sample peaks
    real_unique_peaks_filtered = exclude_background_peaks(unique_peaks_sample, unique_peaks_background,tol=tol)
    print(real_unique_peaks_filtered)

    # Convert to numpy matrix and save
    matrix_maldi = [convert_maldi2numpy(md, target_peaks=real_unique_peaks_filtered, tol=tol, ncpu=ncpu) for md in maldi_datasets]
    [save_numpy_data("/scratch/gpfs/cm7897/output/"+name+".npz",img_data=Z, peaks=real_unique_peaks_filtered) for name,Z in zip(file_names,matrix_maldi)]

    print("Done")
    
    
    
