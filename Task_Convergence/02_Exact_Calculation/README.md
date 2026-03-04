# 02_Exact_Calculation - 精确解计算

## 功能说明

本文件夹包含用于计算热容海平面 (TSLA) 精确解的脚本。精确解直接使用 TEOS-10 非线性状态方程计算密度，不进行任何泰勒展开近似，作为验证近似方法的"标准答案"。

## 文件列表

| 文件名 | 数据源 | 参考态 | 输出 |
|--------|--------|--------|------|
| `Calc_Exact_EN4_Avg.m` | EN4 | 平均态 | `EN4_Formula11_Exact_Avg.mat` |
| `Calc_Exact_EN4_Std.m` | EN4 | 标准态 (S=35, T=0) | `EN4_Formula11_Exact_Std.mat` |
| `Calc_Exact_IAP_Avg.m` | IAP | 平均态 | `IAP_Formula11_Exact_Avg.mat` |
| `Calc_Exact_IAP_Std.m` | IAP | 标准态 | `IAP_Formula11_Exact_Std.mat` |
| `Calc_Exact_Ishii_Avg.m` | Ishii | 平均态 | `Ishii_Formula11_Exact_Avg.mat` |
| `Calc_Exact_Ishii_Std.m` | Ishii | 标准态 | `Ishii_Formula11_Exact_Std.mat` |

## 参考态说明

- **平均态 (Avg)**: 以当地的时间平均 T/S 为展开中心，收敛快，物理意义直观
- **标准态 (Std)**: 全球统一使用 S=35, θ=0 为展开中心，基准统一但需要更高阶展开

## 运行顺序

1. 确保已运行 `01_Data_Preprocessing` 中的平均态计算脚本
2. 可并行运行各数据源的精确解计算

## 依赖

- GSW (TEOS-10) 海洋工具箱
- 平均态文件: `*_Mean_State.mat`

## 注意事项

- 计算精确解需要逐月读取 20 年的全球 3D 数据，计算量巨大
- 单个脚本通常需要 30-60 分钟完成
