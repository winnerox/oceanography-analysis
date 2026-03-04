import os
import numpy as np
import pandas as pd
from scipy.io import loadmat, savemat
import h5py
from tabulate import tabulate

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

def calculate_statistics(data, ds_name, state_cn, var_name, order):
    """计算趋势数据的统计指标"""
    valid_data = data[~np.isnan(data)]
    
    # 基础元数据
    base_info = {
        "数据集": ds_name, "状态": state_cn, "变量": var_name, "阶数": order,
        "类别": f"{ds_name} {var_name}趋势 {order}阶 ({state_cn})"
    }
    
    if len(valid_data) == 0:
        stats = {
            "最小": np.nan, "最大": np.nan, "2.5%": np.nan, "25%": np.nan, 
            "中位数": np.nan, "75%": np.nan, "97.5%": np.nan, "均值": np.nan, "标准差": np.nan
        }
    else:
        stats = {
            "最小": np.min(valid_data), "最大": np.max(valid_data),
            "2.5%": np.percentile(valid_data, 2.5), "25%": np.percentile(valid_data, 25),
            "中位数": np.median(valid_data), "75%": np.percentile(valid_data, 75),
            "97.5%": np.percentile(valid_data, 97.5), "均值": np.mean(valid_data),
            "标准差": np.std(valid_data)
        }
    return {**base_info, **stats}

def print_grouped_table(df, ds_name, state_cn):
    """打印高度优化的连贯表格"""
    print(f"\n\n{'='*20} {ds_name} - {state_cn} 趋势统计 {'='*20}")
    
    # 筛选当前数据集和状态的数据
    sub_df = df[(df['数据集'] == ds_name) & (df['状态'] == state_cn)].copy()
    if sub_df.empty:
        print("（无数据）")
        return

    # 为了终端展示的美观，只挑选最关键的列
    display_cols = ['变量', '阶数', '最小', '最大', '2.5%', '中位数', '97.5%', '均值', '标准差']
    
    # 按照 变量 -> 阶数 排序 (HSLA先显示，然后SSLA, TSLA)
    sub_df['变量'] = pd.Categorical(sub_df['变量'], categories=['HSLA', 'SSLA', 'TSLA', 'Cross'], ordered=True)
    sub_df = sub_df.sort_values(by=['变量', '阶数'])

    # 打印精美表格
    print(tabulate(sub_df[display_cols], headers='keys', tablefmt='rounded_outline', floatfmt=".4f", showindex=False))

def main():
    print("========== 趋势数据高可读性统计表 ==========\n")
    
    trend_dir = r'D:\work\Task_Convergence\Trend_Results'
    datasets = ['EN4', 'IAP', 'Ishii']
    states = {'平均态': 'Average', '标准态': 'StdRef'}

    print('正在加载并计算趋势数据，请稍候...')
    all_stats = []

    for ds_name in datasets:
        for state_cn, state_suffix in states.items():
            try:
                # 构建趋势文件路径
                trend_file = os.path.join(trend_dir, f'{ds_name}_{state_suffix}_Trends.mat')
                
                # 加载趋势数据
                data = load_mat_robust(trend_file)
                
                # 提取趋势数据
                trend_TSLA = data['trend_TSLA']
                trend_SSLA = data['trend_SSLA']
                trend_HSLA = data['trend_HSLA']
                trend_Cross = data['trend_Cross']
                
                # 循环 1 到 8 阶，收集所有数据
                for order in range(1, 9):
                    idx = order - 1
                    all_stats.append(calculate_statistics(trend_HSLA[..., idx], ds_name, state_cn, 'HSLA', order))
                    all_stats.append(calculate_statistics(trend_SSLA[..., idx], ds_name, state_cn, 'SSLA', order))
                    all_stats.append(calculate_statistics(trend_TSLA[..., idx], ds_name, state_cn, 'TSLA', order))
                    all_stats.append(calculate_statistics(trend_Cross[..., idx], ds_name, state_cn, 'Cross', order))
            except FileNotFoundError:
                print(f"警告：未找到文件 {trend_file}")
                pass
            except KeyError as e:
                print(f"警告：文件 {trend_file} 缺少键 {e}")
                pass

    df_all = pd.DataFrame(all_stats)

    # ---------------- 打印极度舒适的连贯表格 ----------------
    for ds_name in ['EN4', 'IAP', 'Ishii']:
        for state_cn in ['平均态', '标准态']:
            print_grouped_table(df_all, ds_name, state_cn)

    # ---------------- 导出至 Excel ----------------
    out_dir = r'D:\work\Task_Convergence'
    os.makedirs(out_dir, exist_ok=True)
    excel_file = os.path.join(out_dir, 'Trend_Statistics_Analysis_Report.xlsx')
    
    if not df_all.empty:
        # 导出 Excel 时，保留所有原始列
        df_all.to_excel(excel_file, index=False)
        print(f'\n✨ 趋势数据统计完毕！完整的 Excel 报告保存在: {excel_file}')

if __name__ == "__main__":
    main()