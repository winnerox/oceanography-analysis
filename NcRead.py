import xarray as xr
import matplotlib.pyplot as plt

# 快速读取
file_path = r"D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc"
ds = xr.open_dataset(file_path)

# 查看数据结构
print(ds)

# 查看变量
print("\n变量列表:", list(ds.data_vars))

# 如果是盐度数据，假设变量名是's'或'salinity'
if 's' in ds:
    # 绘制表层盐度
    ds['s'].isel(lev=0).plot()
    plt.title('表层盐度 (2005年1月)')
    plt.show()
    
    # 查看某个剖面的数据
    profile = ds['s'].isel(lat=50, lon=100)  # 选择特定点
    profile.plot(y='lev')
    plt.gca().invert_yaxis()
    plt.show()

ds.close()