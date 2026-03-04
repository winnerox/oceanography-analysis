import h5py

# 加载 MAT 文件
with h5py.File('d:\\work\\Ishii_mat_data\\Ishii_CrossDetail_Avg.mat', 'r') as f:
    # 打印变量列表
    print('文件中的变量:')
    for key in f.keys():
        print(f'  - {key}')
    
    # 打印关键变量的大小
    print('\n关键变量大小:')
    if 'lon' in f:
        print(f'  lon: {f["lon"].shape}')
    if 'lat' in f:
        print(f'  lat: {f["lat"].shape}')
    if 'time_vec' in f:
        print(f'  time_vec: {f["time_vec"].shape}')
    if 'Cross_T1S1' in f:
        print(f'  Cross_T1S1: {f["Cross_T1S1"].shape}')

print('\n检查完成!')