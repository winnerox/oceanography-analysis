# 04_CrossDetail_Calculation - 详细混合项计算

## 功能说明

本文件夹包含用于计算每一阶混合项细节的脚本，用于详细分析各混合项的衰减规律和物理意义。

## 文件列表

| 文件名 | 数据源 | 参考态 | 阶数范围 |
|--------|--------|--------|----------|
| `Calc_CrossDetail_EN4_Avg.m` | EN4 | 平均态 | 2-5阶 |
| `Calc_CrossDetail_EN4_Std.m` | EN4 | 标准态 | 2-8阶 |
| `Calc_CrossDetail_IAP_Avg.m` | IAP | 平均态 | 2-8阶 |
| `Calc_CrossDetail_IAP_Std.m` | IAP | 标准态 | 2-8阶 |
| `Calc_CrossDetail_Ishii_Avg.m` | Ishii | 平均态 | 2-8阶 |
| `Calc_CrossDetail_Ishii_Std.m` | Ishii | 标准态 | 2-8阶 |

## 混合项说明

混合项是泰勒展开中温度和盐度的交叉项，形如：

- **2阶**: T×S
- **3阶**: T²×S, T×S²
- **4阶**: T³×S, T²×S², T×S³
- ...以此类推

## 用途

- 分析每一阶混合项的空间分布特征
- 验证混合项随阶数增加的衰减规律
- 评估忽略高阶混合项带来的误差

## 运行顺序

1. 确保已运行 `01_Data_Preprocessing` 和 `02_Exact_Calculation`
2. 按需运行各数据源的计算脚本

## 依赖

- GSW (TEOS-10) 海洋工具箱
- TEOS10_HighOrder_Engine (高阶导数计算引擎)
- 平均态文件: `*_Mean_State.mat`
