import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import h5py

# ----------------- 设置中文字体与负号显示 -----------------
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def load_mat_robust(filepath):
    """自适应读取 .mat 文件，兼容旧版 (v7) 和新版 (v7.3 HDF5)"""
    try:
        mat = loadmat(filepath)
        return {k: v for k, v in mat.items() if not k.startswith('__')}
    except NotImplementedError:
        data = {}
        with h5py.File(filepath, 'r') as f:
            for k in f.keys():
                if not k.startswith('#'): 
                    # 注意：HDF5 读取时维度是倒置的，.T 帮你转置回 MATLAB 的顺序
                    data[k] = np.array(f[k]).T
        return data

def plot_spatial_map(data_2d, title, save_path):
    """
    绘制空间分布对比图：左图为原始绝对场，右图为扣除全场平均后的距平场
    """
    # 过滤掉 NaN 值计算统计量
    valid_data = data_2d[~np.isnan(data_2d)]
    spatial_mean = np.mean(valid_data)
    
    # 构造扣除空间均值后的距平场
    data_anomaly = data_2d - spatial_mean

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # --- 左图：原始绝对一阶平均态 ---
    # 使用 RdYlBu_r 颜色盘 (红暖蓝冷)
    im0 = axes[0].imshow(data_2d, cmap='RdYlBu_r', origin='lower', aspect='auto')
    axes[0].set_title(f'{title}\n(原始绝对值, 均值: {spatial_mean:.2f} mm)', fontsize=12)
    fig.colorbar(im0, ax=axes[0], label='TSLA 1阶 (mm)')

    # --- 右图：扣除空间均值后的距平态 ---
    # 以 0 为中心，更能体现数值的相对起伏
    im1 = axes[1].imshow(data_anomaly, cmap='RdBu_r', origin='lower', aspect='auto')
    axes[1].set_title(f'{title}\n(扣除全场平均后的距平场)', fontsize=12)
    fig.colorbar(im1, ax=axes[1], label='TSLA 1阶距平 (mm)')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"✨ 已生成空间分布图: {save_path}")

def main():
    print("========== 开始绘制 1阶平均态空间分布图 ==========\n")
    
    # 1. 定义你的数据文件路径
    files = {
        'EN4': r'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms_1to8_Average.mat',
        'IAP': r'D:\work\IAP_mat_data\IAP_TSLA_Terms_1to8_Average.mat',
        'Ishii': r'D:\work\Ishii_mat_data\Ishii_TSLA_Terms_1to8_Average.mat'
    }
    
    out_dir = r'D:\work\Task_Convergence\Plots'
    os.makedirs(out_dir, exist_ok=True)

    for ds_name, filepath in files.items():
        print(f"正在处理 {ds_name} 数据集...")
        try:
            # 加载数据
            data_dict = load_mat_robust(filepath)
            
            # 提取 TSLA_AllOrders 变量
            if 'TSLA_AllOrders' in data_dict:
                tsla_all = data_dict['TSLA_AllOrders']
                
                # 提取第 1 阶空间数据 (注意：Python 索引从 0 开始)
                # 假设你的数据维度是 (lon, lat, order) 或 (lat, lon, order)
                tsla_1st_order = tsla_all[..., 0] 
                
                # 绘图并保存
                save_path = os.path.join(out_dir, f'{ds_name}_TSLA_1st_Order_Spatial.png')
                plot_spatial_map(tsla_1st_order, f'{ds_name} - TSLA 1阶平均态空间分布', save_path)
            else:
                print(f"  [警告] 在 {ds_name} 中找不到变量 'TSLA_AllOrders'")
                
        except FileNotFoundError:
            print(f"  [错误] 找不到文件: {filepath}")

    print("\n>>> 🎉 所有空间分布图绘制完成！请去 Plots 文件夹查看。 <<<")

if __name__ == "__main__":
    main()