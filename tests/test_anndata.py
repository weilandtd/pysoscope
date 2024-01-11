import numpy as np
import multiprocessing as mp 

from pysoscope.io.anndata import convert_maldi_image_to_anndata, join_anndata
from pysoscope.io.metadata import load_metadata_csv
from pysoscope.io.mat import load_maldi_file_as_dict
from pysoscope.untargeted.sample import subsample_maldi_datasets_naive
from pysoscope.untargeted.pick_peaks import get_unique_peaks
from pysoscope.untargeted.background import get_background_signature

    
if __name__ == "__main__":
    tol = 15e-6  # Peak picking tolerance default 15 ppm
    ncpu = 36
    file_names =["EC2_1","EC2_2",]

    path = "./test_data/"

    metadata = load_metadata_csv("./test_data/metadata.csv",  file_names)

    maldi_files = [path+file+".mat" for file in file_names]

    with mp.Pool(ncpu) as pool:
        maldi_datasets = pool.map(load_maldi_file_as_dict, maldi_files)
    
    # Subsample peaks and get unique peakls
    subsample = subsample_maldi_datasets_naive(maldi_datasets)

        # Get unique sample peaks
    untargeted_peaks = get_unique_peaks(subsample, tol=tol, min_samples=1)

    sample_signature_peaks = [124.007, ]

    background_sig = get_background_signature(maldi_datasets, untargeted_peaks,  sample_signature_peaks,
                                              sig_quantile=0.75,background_quantile=0.25, tol=tol)
    
    maldi_anndatasets = [convert_maldi_image_to_anndata(md, target_peaks=untargeted_peaks, tol=tol, ncpu=ncpu, metadata=metadata[i]) 
                         for i,md in enumerate(maldi_datasets)]
    
    combined_anndata = join_anndata(maldi_anndatasets)

    combined_anndata.write_h5ad('test_anndata.h5ad')