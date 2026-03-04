# 03_Terms_Calculation - 泰勒展开项计算

## 功能说明

本文件夹包含用于计算泰勒展开各项的脚本，包括纯温度项 (Pure T)、纯盐度项 (Pure S) 和混合项 (Cross Terms)。

## 文件列表

### 3.1 纯温度项 (Pure T Terms)

| 文件名 | 数据源 | 参考态 |
|--------|--------|--------|
| `Calc_T_Terms_IAP_Avg.m` | IAP | 平均态 |
| `Calc_T_Terms_IAP_Std.m` | IAP | 标准态 |
| `Calc_T_Terms_Ishii_Avg.m` | Ishii | 平均态 |
| `Calc_T_Terms_Ishii_Std.m` | Ishii | 标准态 |

### 3.2 纯盐度项 (Pure S Terms)

| 文件名 | 数据源 | 参考态 |
|--------|--------|--------|
| `Calc_S_Terms_IAP_Avg.m` | IAP | 平均态 |
| `Calc_S_Terms_IAP_Std.m` | IAP | 标准态 |
| `Calc_S_Terms_Ishii_Avg.m` | Ishii | 平均态 |
| `Calc_S_Terms_Ishii_Std.m` | Ishii | 标准态 |

### 3.3 综合计算

| 文件名 | 功能 | 输出 |
|--------|------|------|
| `Calc_Cross_and_S_Terms_Avg.m` | 计算平均态下的纯S项 + 混合项 (1-3阶) | `EN4_Terms_S_Cross_Avg_1to3.mat` |
| `Calc_Std_Terms_HighOrder.m` | 计算标准态下的高阶项 (1-8阶) | `EN4_Std_Terms_HighOrder_1to8.mat` |

## 泰勒展开公式

TSLA 的泰勒展开可表示为：

```
TSLA ≈ Σ (T项) + Σ (S项) + Σ (混合项)
```

其中混合项形如: T^n × S^m (n+m = 阶数)

## 运行顺序

1. 确保已运行 `01_Data_Preprocessing` 和 `02_Exact_Calculation`
2. 按需运行各数据源和参考态的计算脚本

## 依赖

- GSW (TEOS-10) 海洋工具箱
- TEOS10_HighOrder_Engine (高阶导数计算引擎)
- 平均态文件: `*_Mean_State.mat`
