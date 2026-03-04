import h5py

# 加载 MAT 文件
with h5py.File('d:\\work\\Ishii_mat_data\\Ishii_CrossDetail_Avg.mat', 'r') as f:
    # 打印变量列表和大小
    print('文件中的所有变量及其大小:')
    for key in f.keys():
        print(f'  - {key}: {f[key].shape}')

print('\n检查完成!')