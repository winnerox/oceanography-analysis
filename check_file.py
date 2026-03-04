import os

# 检查文件是否存在
file_path = r"D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc"
print(f"检查文件: {file_path}")
print(f"文件存在: {os.path.exists(file_path)}")

# 检查目录
data_root = r"D:\work\IAP_05_24"
print(f"\n检查目录: {data_root}")
print(f"目录存在: {os.path.exists(data_root)}")

# 列出目录内容
if os.path.exists(data_root):
    print("\n目录内容:")
    for item in os.listdir(data_root):
        print(f"  {item}")