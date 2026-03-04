import scipy.io
import h5py
import os

file_path = 'D:/work/MAT_Data/EN4_TSLA_Terms_1to8_Average.mat'

if not os.path.exists(file_path):
    print(f"Error: File {file_path} not found.")
    exit(1)

try:
    mat = scipy.io.loadmat(file_path)
    print("Loaded with scipy.io")
    for key in mat:
        if not key.startswith('__'):
            val = mat[key]
            print(f"{key}: {val.shape} {val.dtype}")
except NotImplementedError:
    print("scipy.io failed (likely v7.3), trying h5py")
    try:
        with h5py.File(file_path, 'r') as f:
            print("Loaded with h5py")
            for key in f.keys():
                val = f[key]
                print(f"{key}: {val.shape}")
    except Exception as e:
        print(f"h5py failed: {e}")
except Exception as e:
    print(f"scipy.io failed: {e}")
