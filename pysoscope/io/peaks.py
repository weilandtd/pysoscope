import numpy as np

def extract_mz_values(filepath):
    """Extracts m/z values from a file.
        The file is expected to have m/z values in the first column, separated by semicolons.
        Since it is based on exports from SCiLS Lab, initial rows are skipped if they start with a hash (#) (e.g., comments or headers).
    Args:
        filepath (str): Path to the file containing m/z values.
    Returns:
        np.ndarray: Array of m/z values extracted from the file.
        If no valid m/z values are found, an empty array is returned.
    """
    mz_values = []
    
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            mz_value = line.split(";")[0]

            try:
                mz_values.append(float(mz_value))
            except ValueError:
                # Skip lines that don't contain a valid m/z value
                continue

    # Convert list to numpy array
    mz_array = np.array(mz_values)
    return mz_array

def consolidate_peaks(mz_arrays, tol=15e-6, min_samples=2):
    """Consolidates m/z peaks from multiple arrays.
    Peaks are considered the same if they are within a relative tolerance of each other.
    Args:
        mz_arrays (list of np.ndarray): List of arrays containing m/z values.
        tol (float): Relative tolerance for peak consolidation. Defaults to 15e-6.
        min_samples (int): Minimum number of samples required to consider a peak valid. Defaults to 2.
    Returns:
        np.ndarray: Consolidated array of m/z values.
    """
    all_peaks = []
    for peaks in mz_arrays:
        all_peaks.extend(peaks)
    
    all_peaks = np.array(all_peaks)
    all_peaks.sort()
    
    consolidated_peaks = []
    current_peak = all_peaks[0]
    count = 1

    # Iterate through sorted peaks and consolidate
    for peak in all_peaks[1:]:
        # Check if the peak is within the tolerance of the current peak
        if abs(peak - current_peak) <= peak * tol:
            count += 1
        else:
            # If the count of the current peak meets the minimum samples, add it to the consolidated list
            if count >= min_samples:
                consolidated_peaks.append(current_peak)
            current_peak = peak
            count = 1
    
    if count >= min_samples:
        consolidated_peaks.append(current_peak)
    
    return np.array(consolidated_peaks)