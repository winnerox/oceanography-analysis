# 05_Convergence_Analysis - 收敛性分析

## 功能说明

本文件夹包含用于验证泰勒展开收敛性的分析脚本，通过对比精确解与近似解来评估近似方法的准确性。

## 文件列表

| 文件名 | 参考态 | 近似阶数 | 功能 |
|--------|--------|----------|------|
| `Analysis_Convergence.m` | 平均态 | 3阶 | 计算残差分布和收敛性箱线图 |
| `Analysis_Convergence_Std.m` | 标准态 | 8阶 | 验证远离展开中心时的收敛能力 |
| `Analysis_Convergence_Loop_Avg.m` | 平均态 | 循环验证 | 多阶数循环收敛性验证 |
| `Analysis_Convergence_Loop_Std.m` | 标准态 | 循环验证 | 多阶数循环收敛性验证 |

## 分析方法

### 残差计算

```
Residual = Exact - Approx
```

其中 Approx = T(1-n阶) + S(1-n阶) + Cross(2-n阶)

### 输出图表

1. **残差空间分布图**: 展示残差趋势和振幅的全球分布
2. **收敛性箱线图**: 对比 Exact、Approx、Residual 的统计分布

## 核心科学问题

- 泰勒展开保留到第 3 阶是否足够？
- 残差是否已收敛至可忽略的程度？

## 运行顺序

1. 确保已完成 `02_Exact_Calculation` 和 `03_Terms_Calculation`
2. 运行对应参考态的分析脚本

## 依赖

- 精确解数据: `*_Formula11_Exact_*.mat`
- 泰勒展开项数据: `*_Terms_*.mat`
- m_map 工具箱 (地图投影)
