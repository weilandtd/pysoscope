import numpy as np
import anndata as ad


def find_background_pixels(anndata_obj, sample_peaks, percentile, tol=15e-6):
    
    peaks = anndata_obj.var["mz"] # Get the m/z values of the peaks

    # Find the indices of the peaks that are in the sample peaks using tolerance
    sample_peak_indices = []
    for sample_peak in sample_peaks:
        sample_peak_indices.append(np.where(np.abs(peaks - sample_peak) < tol*sample_peak)[0][0])  
    
    # Get the intensities of the sample peaks   
    sample_peak_intensities = anndata_obj.X[:,sample_peak_indices]

    # Find percentiles of the sample peaks
    sample_peak_percentiles = np.percentile(sample_peak_intensities, percentile, axis=0)

    # Find the background pixels as the pixels where all the sample peak intensities are below the percentile
    background_pixels = np.all(sample_peak_intensities < sample_peak_percentiles, axis=1)
    
    return background_pixels


def find_background_signature(anndata_obj, background_pixels, percentile):

    background_signature = np.percentile(anndata_obj.X[background_pixels,:], percentile, axis=0)

    return background_signature


def subtract_background(anndata_obj, sample_peaks, percentile_sample, percentile_background, tol=15e-6):

    modified_anndata_obj = anndata_obj.copy()

    background_pixels = find_background_pixels(anndata_obj, sample_peaks, percentile_sample, tol=tol)

    background_signature = find_background_signature(anndata_obj, background_pixels, percentile_background)

    modified_anndata_obj.X = np.clip(modified_anndata_obj.X - background_signature, 0, np.inf)

    return modified_anndata_obj

