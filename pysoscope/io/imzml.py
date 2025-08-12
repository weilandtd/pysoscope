import xml.etree.ElementTree as ET
import numpy as np
import os
import pandas as pd

def load_imzml_file_as_dict(imzml_filepath):
    """Loads an imzML file and extracts relevant data into a dictionary.
    Args:
        imzml_filepath (str): Path to the imzML file.
    Returns:
        dict: Dictionary containing the following keys:
            - 'data': Dictionary with keys 'id', 'length', 'offset', 'offset2', 'peak_mz', 'peak_sig', 'x', 'y'.
            - 'fname': List containing the filename of the imzML file.
            - 'res': Numpy array containing the resolution.
    Raises:
        FileNotFoundError: If the associated .ibd file is not found.
        ValueError: If the resolution cannot be determined from the metadata and the user does not provide a valid input.
    """
    # Dict template
    maldi_data = {"data": {}, "fname": [], "res": []}

    # Extract filename and related .ibd file
    path, name = os.path.split(imzml_filepath)
    ibd_filepath = os.path.join(path, name.replace('.imzML', '.ibd'))

    if not os.path.exists(ibd_filepath):
        raise FileNotFoundError(f"Associated .ibd file not found: {ibd_filepath}")

    # Parse the imzML XML file
    tree = ET.parse(imzml_filepath)
    root = tree.getroot()

    # Namespace handling (imzML files often use namespaces)
    ns = {'ns': 'http://psi.hupo.org/ms/mzml'}

    # Extract metadata (resolution, etc.)
    scan_settings = root.find('.//ns:scanSettingsList/ns:scanSettings', ns)
    pixel_size_x = scan_settings.find(".//ns:cvParam[@name='pixel size x']", ns)
    pixel_size_y = scan_settings.find(".//ns:cvParam[@name='pixel size y']", ns)

    if pixel_size_x is not None:
        maldi_data['res'] = np.array([[float(pixel_size_x.attrib['value'])]])
    elif pixel_size_y is not None:
        maldi_data['res'] = np.array([[float(pixel_size_y.attrib['value'])]])
    else:
        # Prompt user to supply resolution if not found
        user_res = input("Resolution not found in metadata. Please enter the resolution: ")
        try:
            maldi_data['res'] = np.array([[float(user_res)]])
        except ValueError:
            raise ValueError("Invalid input! Resolution must be a numeric value.")

    maldi_data['fname'] = name

    # Extract spectrum data (x, y positions, and binary data offsets)
    spectra = root.findall('.//ns:spectrum', ns)

    # Read data into the dictionary
    maldi_data['data']['id'] = [np.array([float(i)]) for i in range(1, len(spectra)+1)]
    maldi_data['data']['length'] = []
    maldi_data['data']['offset'] = []
    maldi_data['data']['offset2'] = []
    maldi_data['data']['peak_mz'] = []
    maldi_data['data']['peak_sig'] = []
    maldi_data['data']['x'] = []
    maldi_data['data']['y'] = []

    for spectrum in spectra:
        # Get the x and y positions
        scan = spectrum.find('.//ns:scan', ns)
        position_x = scan.find(".//ns:cvParam[@name='position x']", ns).attrib['value']
        position_y = scan.find(".//ns:cvParam[@name='position y']", ns).attrib['value']
        
        # Get binary data array offsets
        binary_data_arrays = spectrum.findall('.//ns:binaryDataArray', ns)

        # m/z and intensity binary data
        mz_array = binary_data_arrays[0].find(".//ns:cvParam[@name='external offset']", ns).attrib['value']
        mz_length = binary_data_arrays[0].find(".//ns:cvParam[@name='external array length']", ns).attrib['value']
        int_array = binary_data_arrays[1].find(".//ns:cvParam[@name='external offset']", ns).attrib['value']

        # Store the parsed values
        maldi_data['data']['x'].append(float(position_x))
        maldi_data['data']['y'].append(float(position_y))
        maldi_data['data']['length'].append(int(mz_length))
        maldi_data['data']['offset'].append(int(mz_array))
        maldi_data['data']['offset2'].append(int(int_array))

    temp = {"data": {}}
    temp['data']['length'] = np.array(maldi_data['data']['length'])
    temp['data']['offset'] = np.array(maldi_data['data']['offset'])
    temp['data']['offset2'] = np.array(maldi_data['data']['offset2'])

    # Convert to lists of single element numpy arrays
    maldi_data['data']['x'] = [np.array([float(num)]) for num in maldi_data['data']['x']]
    maldi_data['data']['y'] = [np.array([float(num)]) for num in maldi_data['data']['y']]
    maldi_data['data']['length'] = [np.array([float(num)]) for num in maldi_data['data']['length']]
    maldi_data['data']['offset'] = [np.array([float(num)]) for num in maldi_data['data']['offset']]
    maldi_data['data']['offset2'] = [np.array([float(num)]) for num in maldi_data['data']['offset2']]

    # if CSV file with  X,Y coordinates is present, read it
    imzml_name_no_ext = name.replace('.imzML', '')
    csv_filename = None

    # Search for a matching .csv file ignoring case
    for f in os.listdir(path):
        if f.lower().endswith('.csv') and os.path.splitext(f)[0].lower() == imzml_name_no_ext.lower():
            csv_filename = f
            break

    if csv_filename:
        csv_filepath = os.path.join(path, csv_filename)
        # Read the CSV file
        df_temp = pd.read_csv(csv_filepath, skiprows=8, header=0, delimiter=';')
        if 'x' in df_temp.columns and 'y' in df_temp.columns:
            maldi_data['data']['x'] = [np.array([float(x)]) for x in df_temp['x']]
            maldi_data['data']['y'] = [np.array([float(y)]) for y in df_temp['y']]

    # Open and read the ibd file
    with open(ibd_filepath, 'rb') as ibd_file:
        mzml = {'mzml': 'http://psi.hupo.org/ms/mzml'}
        file_type_elem = root.find('.//mzml:cvParam[@accession="IMS:1000031"]', mzml)
        if file_type_elem is not None:
            file_type = file_type_elem.get('name')  # Get the 'name' attribute
            print(f"File type category: {file_type}")
        else:
            print("File type category not found.")

        continuous = file_type=="continuous"

        if continuous:
            raise ValueError("Data type error! This function works only on processed data.")
        else:
            # For processed format, read m/z and intensity arrays for each spectrum
            maldi_data['data']['peak_mz'] = []
            maldi_data['data']['peak_sig'] = []

            for offset, length, offset2 in zip(temp['data']['offset'], temp['data']['length'], temp['data']['offset2']):
                # Read m/z array
                ibd_file.seek(offset)

                bytes_to_read = length * 8  # Convert length to bytes
                mz_data = ibd_file.read(bytes_to_read)
                mz_values = np.float32(np.frombuffer(mz_data, dtype='<f8'))
                maldi_data['data']['peak_mz'].append(mz_values)

                # Read intensity array
                ibd_file.seek(offset2)
                
                bytes_to_read = length * 4  # Convert length to bytes
                intensity_data = ibd_file.read(bytes_to_read)
                intensity_values = np.frombuffer(intensity_data, dtype='<f4')
                maldi_data['data']['peak_sig'].append(intensity_values)

    return maldi_data