# IAP TSLA Terms 数据目录

## 概述
本目录包含IAP（中国科学院大气物理研究所）数据集的泰勒级数展开项计算结果。这些数据用于分析海平面异常（TSLA）的高阶泰勒展开收敛性，与EN4数据集进行对比分析。

## 数据集信息
- **数据源**: IAP（Institute of Atmospheric Physics）海洋温盐再分析数据
- **时间范围**: 2005-2024（20年）
- **空间范围**: 全球网格
- **垂直层**: 可能包含多层深度数据
- **特点**: 中国自主研发的海洋再分析产品，具有较高的空间分辨率

## 文件列表及说明

### 1. 泰勒展开项文件
#### 平均态（Average State）
- `IAP_TSLA_Terms_1to8_Average.mat`: TSLA总展开项（1-8阶）
  - 变量: `TSLA_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 包含TSLA相对于平均态（2005-2024年平均）的1-8阶泰勒展开项

- `IAP_S_Terms_1to8_Average.mat`: 盐度独立项（S-term）
  - 变量: `SSLA_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 盐度变化对TSLA的贡献（独立项）

- `IAP_Cross_Terms_1to8_Average.mat`: 温盐交叉项
  - 变量: `Cross_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 温度和盐度相互作用对TSLA的贡献

#### 标准态（StdRef State）
- `IAP_TSLA_Terms_1to8_StdRef.mat`: TSLA总展开项（1-8阶）
- `IAP_S_Terms_1to8_StdRef.mat`: 盐度独立项
- `IAP_Cross_Terms_1to8_StdRef.mat`: 温盐交叉项

### 2. 公式计算结果
- `IAP_Formula11_Exact_Avg.mat`: 公式11精确值（平均态）
- `IAP_Formula11_Exact_Std.mat`: 公式11精确值（标准态）
  - 变量: 包含`Exact_TSLA`等计算结果

### 3. 平均态文件
- `IAP_Mean_State.mat`: 平均态温盐场
  - 变量: 可能包含`T_mean`, `S_mean`, `P_mean`等

### 4. 其他文件
- `Taylor_Convergence_Test.mat`: 泰勒收敛性测试结果

## 数据格式说明

### MATLAB (.mat) 文件结构
所有MAT文件使用MATLAB v7.3格式（HDF5基础），支持大文件存储。

### 主要变量维度
典型泰勒展开项变量的维度：
```
(time, latitude, longitude, order)
```
- **time**: 时间维度（月/年数据，2005-2024共20年×12月=240个时间点）
- **latitude**: 纬度维度（全球标准网格，可能为1°×1°或更高分辨率）
- **longitude**: 经度维度（全球标准网格）
- **order**: 泰勒展开阶数（1-8阶）

### 数据类型
- 单精度浮点数（float32）或双精度浮点数（float64）
- 无效值用NaN表示

## 物理意义

### 泰勒展开公式
TSLA的泰勒展开表示为：
```
TSLA(t) = TSLA₀ + ∑_{n=1}^8 [1/n! · (∂ⁿTSLA/∂Tⁿ)·ΔTⁿ + 1/n! · (∂ⁿTSLA/∂Sⁿ)·ΔSⁿ + 交叉项]
```

### 各项含义：
1. **TSLA_Terms**: 总展开项 = S_Terms + Cross_Terms + T_Terms
2. **S_Terms**: 仅由盐度变化引起的TSLA变化
3. **Cross_Terms**: 温盐相互作用引起的非线性项

## IAP数据集特点

### 空间分辨率
- 可能为0.5°×0.5°或更高分辨率网格
- 包含边缘海和极区数据

### 数据质量
- 同化多种观测数据（Argo浮标、船测、卫星等）
- 中国近海区域数据质量较高
- 采用先进的数据同化技术

### 时间一致性
- 时间序列连续完整
- 季节循环和长期趋势均被保留

## 使用示例（Python）

```python
import scipy.io as sio
import numpy as np

# 加载IAP TSLA总展开项
data = sio.loadmat('IAP_TSLA_Terms_1to8_Average.mat')
tsla_terms = data['TSLA_AllOrders']  # shape: (time, lat, lon, 8)

# 获取所有网格点的第1阶展开项时间序列
order1_global = tsla_terms[..., 0]  # 第1阶展开项

# 计算全球平均的时间序列
global_mean_ts = np.nanmean(order1_global, axis=(1, 2))
```

## 与EN4数据对比注意事项

1. **网格差异**: IAP和EN4可能使用不同分辨率的网格，需注意插值处理
2. **数据覆盖**: IAP可能在中国近海有更好的数据覆盖
3. **系统偏差**: 不同再分析产品间存在系统偏差，需进行偏差校正

## 相关脚本
- `Analyze_Statistics.py`: 统计分析脚本，计算各阶统计量
- `Compare_Avg_vs_Std.m`: 比较平均态和标准态差异
- `CrossTerm_Trend_Amp_Analysis.m`: 交叉项趋势和振幅分析

## 参考文献
1. IAP数据说明：Cheng et al., 2017
2. TEOS-10: IOC, SCOR and IAPSO, 2010
3. 海洋再分析产品比较研究：Xu et al., 2021