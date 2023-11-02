import numpy as np 
import multiprocessing as mp

from pysoscope.io.mat import load_maldi_file_as_dict as load_maldi_file
from pysoscope.untargeted.sample import subsample_maldi_file
from pysoscope.untargeted.pick_peaks import get_unique_peaks
from pysoscope.untargeted.convert import convert_maldi2numpy, save_numpy_data


if __name__ == "__main__":
    # Parameters and files 

    tol = 15e-6  # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names = ["R1_9AA", "R3_9AA", "R5_9AA", "R7_9AA",
                  "R2_1_9AA","R6_1_9AA","R7_1_9AA", "R2_9AA", "R8_9AA" ,"R4_9AA",]
    path = "./cm_data/"

    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file, maldi_files)
    
    # Subsample peaks and get unique peakls
    subsample = np.concatenate([subsample_maldi_file(md,N=500) for md in maldi_datasets])

    # Get unique sample peaks
    real_unique_peaks_filtered = get_unique_peaks(subsample, tol=tol)
    print(real_unique_peaks_filtered)
    
    # Convert to numpy matrix and save
    matrix_maldi = [convert_maldi2numpy(md, target_peaks=real_unique_peaks_filtered, tol=tol, ncpu=ncpu) for md in maldi_datasets]
    [save_numpy_data("./cm_data/"+name+".npz",img_data=Z, peaks=real_unique_peaks_filtered) for name,Z in zip(file_names,matrix_maldi)]

    print("Done")