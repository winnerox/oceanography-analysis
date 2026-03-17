# 04_plot（绘图与结果分析）

本目录包含所有绘图与结果分析脚本，基于 03_calculate_trend_amplitude 的输出生成空间分布、箱线图、集合平均与对比诊断。

输入
- D:\work\Task_Convergence\Trend_Results\*
- D:\work\Task_Convergence\Amplitude_Results\*

输出
- 默认保存到：D:\work\Figures

脚本说明
- D:\work\Final_version\04_plot\Ensemble_Analysis_AllInOne.m
  - 功能：多数据集集合平均与置信区间分析
  - 支持：
    - 空间插值到统一 1°×1° 网格
    - 集合均值与 90% 置信区间检验
    - TSLA/HSLA/SSLA 目标可切换

- D:\work\Final_version\04_plot\Plot_Trend_Combined.m
  - 功能：趋势空间分布 + 箱线图
  - 特点：
    - 根据显著性打点遮罩
    - 自适应颜色范围
    - 支持 SSLA/TSLA/HSLA

- D:\work\Final_version\04_plot\Plot_Amplitude_Combined.m
  - 功能：年周期振幅空间分布 + 箱线图

- D:\work\Final_version\04_plot\Plot_CrossDetail_Trend.m
  - 功能：交叉项趋势空间分布

- D:\work\Final_version\04_plot\Plot_Cumulative_Sum.m
  - 功能：累积和诊断与趋势分解展示

- D:\work\Final_version\04_plot\Plot_Exact_SSLA_Difference.m
  - 功能：不同数据集 / 状态下 SSLA 精确项对比

- D:\work\Final_version\04_plot\Plot_State_Difference.m
  - 功能：Average 与 StdRef 状态差异对比

- D:\work\Final_version\04_plot\Plot_Steric_Budget.m
  - 功能：Steric budget 分解与综合图

运行建议
- 先完成 02_calculate 与 03_calculate_trend_amplitude
- 再运行此目录脚本生成图件

依赖
- MATLAB
- m_map 工具箱（用于地图投影与海岸线）
- Trend_Results / Amplitude_Results 产物
