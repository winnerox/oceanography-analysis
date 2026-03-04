import h5py
import numpy as np
import os
import sys

def check_h5_file(filepath):
    print(f'\n{"="*60}')
    print(f'Checking: {filepath}')
    print(f'{"="*60}')
    
    try:
        print(f'File exists: {os.path.exists(filepath)}')
        print(f'File size: {os.path.getsize(filepath)/1024/1024:.2f} MB')
        
        with h5py.File(filepath, 'r') as f:
            print(f'\nFile structure:')
            
            def print_structure(name, obj):
                if isinstance(obj, h5py.Group):
                    print(f'  Group: {name}')
                elif isinstance(obj, h5py.Dataset):
                    print(f'  Dataset: {name}')
                    print(f'    shape: {obj.shape}')
                    print(f'    dtype: {obj.dtype}')
                    print(f'    size: {obj.size}')
                    
                    # 读取部分数据进行分析
                    if obj.size > 0:
                        try:
                            # 对于大文件，只读取一小部分
                            if obj.size > 1000000:
                                sample_size = min(1000, obj.size)
                                if len(obj.shape) == 1:
                                    sample = obj[:sample_size]
                                else:
                                    # 取第一个元素
                                    sample = obj[0]
                            else:
                                sample = obj[:]
