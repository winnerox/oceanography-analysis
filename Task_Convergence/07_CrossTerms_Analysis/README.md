# 07_CrossTerms_Analysis - 混合项专项分析

## 功能说明

本文件夹包含用于深入分析混合项 (Cross Terms) 特征的脚本，包括量级对比和空间分布分析。

## 文件列表

| 文件名 | 功能 |
|--------|------|
| `Analysis_CrossTerms.m` | 对比平均态 vs 标准态的混合项总量级 |
| `Analysis_CrossDetail_Avg.m` | 平均态下每个混合项的空间分布 |
| `Analysis_CrossDetail_Std.m` | 标准态下每个混合项的空间分布 |

## 混合项类型

| 阶数 | 混合项 | 物理意义 |
|------|--------|----------|
| 2阶 | T×S | 温盐基本耦合 |
| 3阶 | T²×S, T×S² | 温盐二阶耦合 |
| 4阶 | T³×S, T²×S², T×S³ | 温盐三阶耦合 |
| ... | ... | ... |

## 分析方法

### 量级对比

计算全球 RMS 量级：

```
Magnitude = sqrt(mean(data²))
```

### 空间分布

绘制各混合项的趋势和振幅全球分布图

## 输出图表

1. 混合项总量级柱状对比图
2. 各混合项空间分布图
3. 平均态 vs 标准态对比图

## 运行顺序

1. 确保已完成 `03_Terms_Calculation` 和 `04_CrossDetail_Calculation`
2. 运行分析脚本

## 依赖

- 混合项数据: `*_Terms_S_Cross_*.mat`, `*_CrossDetail_*.mat`
- m_map 工具箱 (地图投影)
