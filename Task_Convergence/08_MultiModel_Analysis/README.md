# 08_MultiModel_Analysis - 多模式验证分析

## 功能说明

本文件夹包含用于跨数据集验证结果稳健性的脚本，对比 EN4、IAP、Ishii 三套不同数据产品的结果。

## 文件列表

| 文件名 | 参考态 | 功能 |
|--------|--------|------|
| `Analysis_MultiModel_Significance_Avg.m` | 平均态 | 多模式显著性分析 |
| `Analysis_MultiModel_Significance_Std.m` | 标准态 | 多模式显著性分析 |

## 数据产品说明

| 数据集 | 来源 | 特点 |
|--------|------|------|
| EN4 | UK Met Office | 经典客观分析数据集 |
| IAP | 中科院大气物理研究所 | 中国自主数据集 |
| Ishii | JMA (日本气象厅) | 日本海洋数据集 |

## 分析方法

### 集合平均

```
Ensemble_Mean = mean(EN4, IAP, Ishii)
```

### 信噪比 (SNR)

```
SNR = Ensemble_Mean / Ensemble_Std
```

SNR > 1 表示信号显著

## 输出图表

1. 集合平均趋势图
2. 信噪比打点图 (显著性标记)
3. 多模式对比箱线图

## 科学意义

- 验证结果对不同数据源的敏感性
- 评估结论的普适性
- 识别数据产品间的差异

## 运行顺序

1. 确保已完成 EN4、IAP、Ishii 三套数据的计算
2. 运行多模式分析脚本

## 依赖

- 各数据源的精确解和泰勒展开项数据
- m_map 工具箱 (地图投影)
