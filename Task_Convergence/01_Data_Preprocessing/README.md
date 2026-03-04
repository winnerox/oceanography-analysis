# 01_Data_Preprocessing - 数据预处理

## 功能说明

本文件夹包含用于计算海洋数据平均态的脚本，为后续的泰勒展开分析提供基础数据。

## 文件列表

| 文件名 | 功能 | 输出 |
|--------|------|------|
| `Calc_Mean_EN4.m` | 计算 EN4 数据的平均态场 (T/S) | `EN4_Mean_State.mat` |
| `Calc_Mean_IAP.m` | 计算 IAP 数据的平均态场 (T/S) | `IAP_Mean_State.mat` |
| `Calc_Mean_Ishii.m` | 计算 Ishii 数据的平均态场 (T/S) | `Ishii_Mean_State.mat` |

## 运行顺序

1. 先运行 `Calc_Mean_EN4.m`
2. 再运行 `Calc_Mean_IAP.m`
3. 最后运行 `Calc_Mean_Ishii.m`

## 依赖

- 输入数据路径: `D:\work\IAP_05_24`, `D:\work\Ishii_05_24`
- 输出路径: `D:\work\IAP_TSLA_Terms`, `D:\work\Ishii_TSLA_Terms`

## 注意事项

- 这些脚本需要读取 20 年的月度数据 (2005-2024)，计算量较大
- 输出的平均态文件是后续所有计算的基础，必须首先运行
