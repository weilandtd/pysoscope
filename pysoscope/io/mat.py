import h5py
import numpy as np


def load_maldi_file_as_dict(filepath,):

    # Dict template
    maldi_data = {"data":{} ,"fname":[], "res":[]}
    
    f = h5py.File(filepath)
    
    # Load meta data
    msi = f['msi']
    data = msi['data']
    maldi_data['fname'] = np.array(msi['fname'])
    maldi_data['res'] = np.array(msi['res'])

    # Read the data
    for k, v in data.items():
        if not ( k in ["peak_mz"], ["peak_sig"]):
            maldi_data["data"][k] = [np.array(f[vi[0]])[0][0] for vi in v ]
        else:
            maldi_data["data"][k] = [np.array(f[vi[0]])[0] for vi in v ]
 
    # Convert X/Z in maldi data 
    _x = maldi_data["data"]["x"]
    _y = maldi_data["data"]["y"]

    min_x = np.min(_x)
    min_y = np.min(_y)

    x = ( _x - min_x) + 1
    y = ( _y - min_y) + 1

    print(x[:,0])
    maldi_data["data"]["x0"] = x[:,0]
    maldi_data["data"]["y0"] = y[:,0]
    
    return maldi_data