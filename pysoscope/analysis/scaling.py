import numpy as np 

def scale_no_zeros(anndata_obj):
    """
    Scales the anndata object by subtracting the mean and dividing by the standard deviation of the peak associated data,
    but excludes the pixels where the peak associated data is zero from scaling, while including them in the final output.
    """

    modified_anndata_obj = anndata_obj.copy()

    data = anndata_obj.X

    # Exclude zeros from scaling
    non_zero_data = data.copy()
    non_zero_data[non_zero_data == 0] = np.nan

    # Calculate mean and standard deviation of non-zero data
    mean = np.nanmean(non_zero_data, axis=0)
    std = np.nanstd(non_zero_data, axis=0)

    # Scale the data by subtracting mean and dividing by standard deviation
    scaled_data = (data - mean) / std

    # Replace the scaled values of zeros with zeros
    scaled_data[data == 0] = 0

    modified_anndata_obj.X = scaled_data

    return modified_anndata_obj

    