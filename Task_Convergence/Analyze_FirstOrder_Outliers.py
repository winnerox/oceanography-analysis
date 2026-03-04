import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import h5py

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def load_mat_robust(filepath):
    try:
        mat = loadmat(filepath)
        return {k: v for k, v in mat.items() if not k.startswith('__')}
    except NotImplementedError:
        data = {}
        with h5py.File(filepath, 'r') as f:
            for k in f.keys():
                if not k.startswith('#'):
                    data[k] = np.array(f[k]).T
        return data

def plot_spatial_map(data_2d, title, save_path):
    # data_2d shape: (Lat, Lon)
    valid_data = data_2d[~np.isnan(data_2d)]
    spatial_mean = np.nanmean(valid_data) if len(valid_data) > 0 else np.nan
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    im0 = ax.imshow(data_2d, cmap='RdYlBu_r', origin='lower', aspect='auto')
    ax.set_title(f'{title}\n(空间均值: {spatial_mean:.2f} mm)', fontsize=14)
    fig.colorbar(im0, ax=ax, label='mm')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    
def plot_boxplot(data_dict, title, save_path):
    labels = [k for k in data_dict.keys() if len(data_dict[k]) > 0]
    data = [data_dict[k] for k in labels]
    
    if not data:
        print(f"      [跳过] {title} 无有效数据")
        return
        
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.boxplot(data, positions=range(1, len(labels)+1), patch_artist=True, showfliers=True, 
               boxprops=dict(facecolor='lightblue', color='k'),
               flierprops=dict(marker='o', color='r', markersize=2, alpha=0.5))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('1阶值 (mm)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

def main():
    print("========== 开始绘制 1阶值分布与箱线图 ==========\n")
    base_dir = r'D:\work'
    dir_map = {
        'EN4': os.path.join(base_dir, 'EN4_TSLA_Terms'),
        'IAP': os.path.join(base_dir, 'IAP_mat_data'),
        'Ishii': os.path.join(base_dir, 'Ishii_mat_data')
    }
    
    out_dir = os.path.join(base_dir, 'Task_Convergence', 'Plots', 'FirstOrder_Outliers')
    os.makedirs(out_dir, exist_ok=True)
    
    models = ['EN4', 'IAP', 'Ishii']
    states = {'Average': '平均态', 'StdRef': '标准态'}
    vars_names = ['HSLA', 'SSLA', 'TSLA']
    
    for state, state_cn in states.items():
        print(f"正在处理 {state_cn}...")
        
        # 准备收集用于箱线图的数据
        var_data = {v: {m: [] for m in models} for v in vars_names}
        
        for ds_name in models:
            print(f"  加载 {ds_name}...")
            d_dir = dir_map[ds_name]
            tsla_file = os.path.join(d_dir, f'{ds_name}_TSLA_Terms_1to8_{state}.mat')
            ssla_file = os.path.join(d_dir, f'{ds_name}_S_Terms_1to8_{state}.mat')
            cross_file = os.path.join(d_dir, f'{ds_name}_Cross_Terms_1to8_{state}.mat')
            
            try:
                # 兼容部分文件可能缺失的情况
                if not (os.path.exists(tsla_file) and os.path.exists(ssla_file) and os.path.exists(cross_file)):
                    print(f"    [跳过] {ds_name} 文件不全")
                    continue
                    
                tsla_all = load_mat_robust(tsla_file)['TSLA_AllOrders'][..., 0] # shape (Lon, Lat, Time)
                ssla_all = load_mat_robust(ssla_file)['SSLA_AllOrders'][..., 0]
                cross_all = load_mat_robust(cross_file)['Cross_AllOrders'][..., 0]
                hsla_all = tsla_all + ssla_all + cross_all
                
                dataset_vars = {'HSLA': hsla_all, 'SSLA': ssla_all, 'TSLA': tsla_all}
                
                for v_name, v_data in dataset_vars.items():
                    # 绘制空间时间平均图
                    v_mean = np.nanmean(v_data, axis=2) # shape (Lon, Lat)
                    # 转置为 (Lat, Lon) for imshow origin='lower'
                    v_mean = v_mean.T
                    
                    save_path = os.path.join(out_dir, f'{ds_name}_{v_name}_1stOrder_{state}_Spatial.png')
                    plot_spatial_map(v_mean, f'{ds_name} {v_name} 1阶各点时间平均 ({state_cn})', save_path)
                    
                    # 收集用于箱线图的数据 (去除 NaN 并展平)
                    valid_data = v_data[~np.isnan(v_data)].flatten()
                    # 为了加速绘图，采样一定数量的点
                    if len(valid_data) > 100000:
                        valid_data = np.random.choice(valid_data, 100000, replace=False)
                    var_data[v_name][ds_name] = valid_data
                    
            except Exception as e:
                print(f"    [错误] 处理 {ds_name} 时报错: {e}")
                
        # 绘制该状态下每个变量的三家数据集对比箱线图
        for v_name in vars_names:
            print(f"  绘制 {v_name} 箱线图...")
            save_path = os.path.join(out_dir, f'{v_name}_1stOrder_{state}_Boxplot.png')
            plot_boxplot(var_data[v_name], f'{v_name} 1阶值分布对比 ({state_cn})', save_path)

    print("\n🎉 分析完成！所有图表均保存在:", out_dir)

if __name__ == '__main__':
    main()
