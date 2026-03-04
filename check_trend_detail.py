import scipy.io as sio
import numpy as np
import sys

try:
    data = sio.loadmat('d:\\work\\Trend_Result_Final.mat')
    print('Keys:', [k for k in data.keys() if not k.startswith('__')])
    
    total_info = data['total_info_Trend'][0, 0]
    print('\nStructure fields:')
    for field in total_info.dtype.names:
        field_data = total_info[field][0, 0] if total_info[field].size > 0 else total_info[field]
        print(f'  {field}: shape={field_data.shape}, dtype={field_data.dtype}')
        
        # 如果是数值数组，显示一些统计信息
        if hasattr(field_data, 'dtype') and np.issubdtype(field_data.dtype, np.number):
            if field_data.size > 0:
                flat_data = field_data.flatten()
                valid_data = flat_data[~np.isnan(flat_data)]
                if len(valid_data) > 0:
                    print(f'    min={valid_data.min():.6f}, max={valid_data.max():.6f}, mean={valid_data.mean():.6f}')
    
    # 特别检查Trend_mean和Trend_sig
    print('\nDetailed Trend_mean info:')
    trend_mean = total_info['Trend_mean'][0, 0]
    print(f'  Trend_mean shape: {trend_mean.shape}')
    if trend_mean.size > 0:
        print(f'  Trend_mean sample values (first 5): {trend_mean.flatten()[:5]}')
    
    print('\nDetailed Trend_sig info:')
    trend_sig = total_info['Trend_sig'][0, 0]
    print(f'  Trend_sig shape: {trend_sig.shape}')
    if trend_sig.size > 0:
        print(f'  Trend_sig sample values (first 5): {trend_sig.flatten()[:5]}')
        
except Exception as e:
    print(f'Error: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)