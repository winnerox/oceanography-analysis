import scipy.io

# 加载 MAT 文件
mat = scipy.io.loadmat('d:\\work\\Ishii_mat_data\\Ishii_CrossDetail_Avg.mat')

# 打印变量列表
print('文件中的变量:')
for key in mat.keys():
    if not key.startswith('__'):
        print(f'  - {key}')

# 打印关键变量的大小
print('\n关键变量大小:')
if 'lon' in mat:
    print(f'  lon: {mat["lon"].shape}')
if 'lat' in mat:
    print(f'  lat: {mat["lat"].shape}')
if 'time_vec' in mat:
    print(f'  time_vec: {mat["time_vec"].shape}')
if 'Cross_T1S1' in mat:
    print(f'  Cross_T1S1: {mat["Cross_T1S1"].shape}')

print('\n检查完成!')