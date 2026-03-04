import scipy.io as sio
import numpy as np
import sys

try:
    data = sio.loadmat('d:\\work\\Trend_Result_Final.mat')
    print('Keys:', [k for k in data.keys() if not k.startswith('__')])
    for k in data.keys():
        if not k.startswith('__'):
            v = data[k]
            print(f'{k}: shape={v.shape}, dtype={v.dtype}')
            if hasattr(v, 'dtype'):
                print(f'   min={np.nanmin(v) if v.size > 0 else "N/A"}, max={np.nanmax(v) if v.size > 0 else "N/A"}')
except Exception as e:
    print(f'Error: {e}')
    sys.exit(1)