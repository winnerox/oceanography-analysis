# 02_calculate（数据集主计算）

本目录包含三套数据集（EN4、IAP、Ishii）的主计算脚本。
核心目标是：在两种参考状态（Average 与 StdRef）下，计算海水密度扰动导致的海面高度变化分解，并输出可用于趋势与绘图的时空矩阵。

总体输出（按数据集分类）
- D:\work\EN4_mat_data\*
- D:\work\IAP_mat_data\*
- D:\work\Ishii_mat_data\*

输出内容包含
- SSLA（Exact）：基于真实密度计算的精确总比容项
- TSLA（Thermosteric）：温度项的泰勒展开 1–8 阶
- HSLA（Halosteric）：盐度项的泰勒展开 1–8 阶
- Cross（混合项）：温度×盐度交叉高阶项
- CrossDetail（细分交叉项）：Cross_TkSj 的逐项结果
- Time_Axis / lon / lat

脚本说明
- D:\work\Final_version\02_calculate\EN4_4in1.m
  - 输入数据目录：D:\work\EN4_analyses_c13_last20years
  - 计算周期：2005–2024（按月）
  - 过程：
    - 构建 0–2000 m 深度积分权重（dz）
    - 计算多年平均状态（Mean）并缓存
    - 调用 TEOS10_General_Engine 计算高阶导数
    - 计算 Average 与 StdRef 两套结果
  - 输出目录：D:\work\EN4_mat_data

- D:\work\Final_version\02_calculate\IAP_4in1.m
  - 输入数据目录：D:\work\IAP_05_24\TEMP 与 D:\work\IAP_05_24\SALT
  - 计算周期：2005–2024（按月）
  - 特点：
    - 盐度经向插值（若温盐经度网格不一致）
    - depth_bnds 缺失时自动重建 dz

- D:\work\Final_version\02_calculate\Ishii_4in1.m
  - 输入数据目录：D:\work\Ishii_05_24\Temperature 与 D:\work\Ishii_05_24\Salinity
  - 计算周期：2005–2024（按年文件，内部分月）
  - 特点：
    - 按年读取后逐月切片处理
    - depth_bnds 缺失时自动重建 dz

计算核心
- 使用 TEOS-10（GSW）函数计算：
  - Absolute Salinity / Conservative Temperature
  - 海水密度与重力加速度
- 通过泰勒展开计算温度项、盐度项、交叉项贡献
- 将 0–2000m 垂向积分转化为海面高度变化（mm）

输出命名规则（示例）
- EN4_SSLA_Exact_Avg.mat
- EN4_TSLA_Terms_1to8_StdRef.mat
- EN4_CrossDetail_Std.mat

运行依赖
- MATLAB
- GSW (TEOS-10) 工具箱
- 原始 NetCDF 数据
