import scipy.io as sio
import numpy as np
import sys
import os

def check_mat_file(filepath):
    print(f'\n{"="*60}')
    print(f'Checking: {filepath}')
    print(f'{"="*60}')
    
    try:
        data = sio.loadmat(filepath)
        print(f'File exists: {os.path.exists(filepath)}')
        print(f'File size: {os.path.getsize(filepath)/1024/1024:.2f} MB')
        
        keys = [k for k in data.keys() if not k.startswith('__')]
        print(f'Keys: {keys}')
        
        for k in keys:
            v = data[k]
            print(f'\n  {k}:')
            print(f'    shape: {v.shape}')
            print(f'    dtype: {v.dtype}')
            print(f'    size: {v.size}')
            
            if hasattr(v, 'dtype') and np.issubdtype(v.dtype, np.number):
                if v.size > 0:
                    flat_data = v.flatten()
                    valid_data = flat_data[~np.isnan(flat_data)]
                    if len(valid_data) > 0:
                        print(f'    valid values: {len(valid_data)}/{v.size}')
                        print(f'    min: {valid_data.min():.6e}')
                        print(f'    max: {valid_data.max():.6e}')
                        print(f'    mean: {valid_data.mean():.6e}')
                        print(f'    std: {valid_data.std():.6e}')
                    
    except Exception as e:
        print(f'Error loading file: {e}')

# 检查EN4的TSLA文件
check_mat_file('d:\\work\\EN4_TSLA_Terms\\EN4_TSLA_Terms_1to8_Average.mat')
check_mat_file('d:\\work\\EN4_TSLA_Terms\\EN4_TSLA_Terms_1to8_StdRef.mat')

# 检查IAP的TSLA文件
check_mat_file('d:\\work\\IAP_TSLA_Terms\\IAP_TSLA_Terms_1to8_Average.mat')
check_mat_file('d:\\work\\IAP_TSLA_Terms\\IAP_TSLA_Terms_1to8_StdRef.mat')

# 检查Ishii的TSLA文件（在Ishii_mat_data目录）
check_mat_file('d:\\work\\Ishii_mat_data\\Ishii_TSLA_Terms_1to8_Average.mat')
check_mat_file('d:\\work\\Ishii_mat_data\\Ishii_TSLA_Terms_1to8_StdRef.mat')