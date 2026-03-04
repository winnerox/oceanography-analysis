# EN4 TSLA Terms 数据目录

## 概述
本目录包含EN4数据集的泰勒级数展开项（Taylor Series Expansion Terms）计算结果。这些数据用于分析海平面异常（TSLA）的高阶泰勒展开收敛性。

## 数据集信息
- **数据源**: EN4 (EN.4.2.2) 海洋温盐再分析数据
- **时间范围**: 2005-2024（20年）
- **空间范围**: 全球网格
- **垂直层**: 可能包含多层深度数据

## 文件列表及说明

### 1. 泰勒展开项文件
#### 平均态（Average State）
- `EN4_TSLA_Terms_1to8_Average.mat`: TSLA总展开项（1-8阶）
  - 变量: `TSLA_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 包含TSLA相对于平均态（2005-2024年平均）的1-8阶泰勒展开项

- `EN4_S_Terms_1to8_Average.mat`: 盐度独立项（S-term）
  - 变量: `SSLA_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 盐度变化对TSLA的贡献（独立项）

- `EN4_Cross_Terms_1to8_Average.mat`: 温盐交叉项
  - 变量: `Cross_AllOrders` - 形状: `(time, lat, lon, 8)`
  - 描述: 温度和盐度相互作用对TSLA的贡献

#### 标准态（StdRef State）
- `EN4_TSLA_Terms_1to8_StdRef.mat`: TSLA总展开项（1-8阶）
- `EN4_S_Terms_1to8_StdRef.mat`: 盐度独立项
- `EN4_Cross_Terms_1to8_StdRef.mat`: 温盐交叉项

### 2. 公式计算结果
- `EN4_Formula11_Exact_Avg.mat`: 公式11精确值（平均态）
- `EN4_Formula11_Exact_Std.mat`: 公式11精确值（标准态）
  - 变量: 包含`Exact_TSLA`等计算结果

### 3. 平均态文件
- `EN4_Mean_State.mat`: 平均态温盐场
  - 变量: 可能包含`T_mean`, `S_mean`, `P_mean`等

### 4. 其他文件
- `TEOS10_Engine_Cache.mat`: TEOS-10状态方程计算缓存
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
- **latitude**: 纬度维度（全球标准网格，如1°分辨率）
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
1. **TSLA_Terms**: 总展开项 = S_Terms + Cross_Terms + T_Terms（T_Terms可能已包含在TSLA_Terms中）
2. **S_Terms**: 仅由盐度变化引起的TSLA变化
3. **Cross_Terms**: 温盐相互作用引起的非线性项

## 使用示例（Python）

```python
import scipy.io as sio
import numpy as np

# 加载TSLA总展开项
data = sio.loadmat('EN4_TSLA_Terms_1to8_Average.mat')
tsla_terms = data['TSLA_AllOrders']  # shape: (time, lat, lon, 8)

# 获取第5阶展开项
order_5 = tsla_terms[..., 4]  # 索引4对应第5阶（0-based）

# 计算空间平均的时间序列
spatial_mean = np.nanmean(order_5, axis=(1, 2))
```

## 相关脚本
- `Analyze_Statistics.py`: 统计分析脚本，计算各阶统计量
- `Test_Taylor_Convergence.m`: 测试泰勒级数收敛性

## 注意事项
1. 数据已进行质量控制，无效值标记为NaN
2. 空间网格可能为规则经纬度网格
3. 时间维度为连续月数据
4. 单位为：TSLA (m), 温度 (°C), 盐度 (PSU)

## 参考文献
1. TEOS-10: IOC, SCOR and IAPSO, 2010
2. EN4数据说明：Good et al., 2013