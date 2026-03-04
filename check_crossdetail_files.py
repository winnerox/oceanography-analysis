import h5py
import numpy as np

# 检查文件结构
def check_file_structure(file_path):
    print(f"\n=== 检查文件: {file_path} ===")
    try:
        with h5py.File(file_path, 'r') as f:
            print("文件结构:")
            def print_structure(name, obj):
                if isinstance(obj, h5py.Group):
                    print(f"  组: {name}")
                elif isinstance(obj, h5py.Dataset):
                    print(f"  数据集: {name} - 形状: {obj.shape}, 类型: {obj.dtype}")
            f.visititems(print_structure)
    except Exception as e:
        print(f"错误: {e}")

# 检查EN4文件
print("=== EN4 CrossDetail文件 ===")
en4_avg_file = "d:/work/EN4_TSLA_Terms/EN4_CrossDetail_Avg.mat"
en4_std_file = "d:/work/EN4_TSLA_Terms/EN4_CrossDetail_Std.mat"

check_file_structure(en4_avg_file)
check_file_structure(en4_std_file)

# 检查Ishii文件
print("\n=== Ishii CrossDetail文件 ===")
ishii_avg_file = "d:/work/Ishii_mat_data/Ishii_CrossDetail_Avg.mat"
ishii_std_file = "d:/work/Ishii_mat_data/Ishii_CrossDetail_Std.mat"

check_file_structure(ishii_avg_file)
check_file_structure(ishii_std_file)
