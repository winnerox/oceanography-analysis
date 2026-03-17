# 01_base（基础工具）

本目录提供整个流程的核心基础工具，供后续计算、趋势/振幅估计与绘图脚本调用。

主要功能
- 时间序列谐波回归：估计线性趋势、年周期与半年周期的幅相与显著性。
- TEOS-10 高阶导数引擎：计算海水热力学性质的高阶纯导数与混合导数。

文件说明
- D:\work\Final_version\01_base\gmt_harmonic_new.m
  - 功能：对时间序列进行最小二乘谐波拟合，模型包含：
    - 常数项（偏置）
    - 线性趋势项
    - 年周期（cos/sin）
    - 半年周期（cos/sin）
  - 输入：
    - t：时间序列（1×N 或 N×1）
    - t1：用于插值/预测的时间序列（可为空）
    - grid_data：2D 或 3D 网格时间序列（lat×lon×time）
    - grid_data_std（可选）：对应标准差
  - 输出：
    - 年/半年振幅与相位（Amplitude1/2, Phase1/2）
    - 线性趋势、标准差、显著性（Trend, Trend_std, Trend_sig, Trend_p, Trend_t）
    - 预测趋势线与拟合序列（Trend_line, Interp）
    - 残差（Resid）
  - 典型用途：为趋势/振幅产品提供逐格点的趋势与显著性判断。

- D:\work\Final_version\01_base\TEOS10_General_Engine.m
  - 功能：TEOS-10 高阶导数引擎，支持任意阶（默认 1–8）纯导数与混合导数的矢量化计算。
  - 特性：
    - 通过系数表与二项式系数加速
    - 可直接输出一整套导数网格
  - 典型用途：在各数据集主计算脚本中快速计算 TSLA/HSLA/Cross 的高阶泰勒项。

在流程中的位置
- 该目录是所有计算与分析脚本的“底座”，不会单独产出结果文件。
- 其函数被以下目录调用：
  - D:\work\Final_version\02_calculate
  - D:\work\Final_version\03_calculate_trend_amplitude

运行依赖
- MATLAB
- GSW (TEOS-10) 工具箱
- 各数据集对应的原始 NetCDF 文件
