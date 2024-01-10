#
def find_duplicates_with_extra_neutron(mz_values, mean_intensities, max_intensity_percent=0.1, tol=15e-6):
    """Find duplicates with N extra neutron."""
    
    duplicates = []
    for i in range(len(mz_values)):
        for j in range(i+1, len(mz_values)):
            if abs(mz_values[i] - mz_values[j]) < 1.003 * tol \
               and abs(mean_intensities[i] - mean_intensities[j]) <= max_intensity_percent * mean_intensities[i]:
                duplicates.append((mz_values[i], mz_values[j]))

    return duplicates
