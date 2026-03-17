# 03_calculate_trend_amplitude（趋势与振幅计算）

本目录负责将 02_calculate 产生的时空序列，转化为可用于图件的趋势与年周期振幅产品。

主要目标
- 对 SSLA / TSLA / HSLA 以及 Cross 逐项进行谐波回归
- 产出：趋势、显著性、年周期振幅
- 输出为 3D 堆叠矩阵，直接服务绘图脚本

输入
- D:\work\EN4_mat_data\*
- D:\work\IAP_mat_data\*
- D:\work\Ishii_mat_data\*

输出
- D:\work\Task_Convergence\Trend_Results\*
- D:\work\Task_Convergence\Amplitude_Results\*

脚本说明
- D:\work\Final_version\03_calculate_trend_amplitude\Build_Trend_Amplitude_Data.m
  - 处理 SSLA 精确项 + TSLA/HSLA 的 1–8 阶
  - 使用 gmt_harmonic_new 进行逐格点回归
  - 输出：
    - trend_SSLA, sig_SSLA, amp_SSLA
    - trend_TSLA, sig_TSLA, amp_TSLA
    - trend_HSLA, sig_HSLA, amp_HSLA
  - 保存到 Trend_Results / Amplitude_Results

- D:\work\Final_version\03_calculate_trend_amplitude\Build_Cross_Trend_Amplitude.m
  - 处理 28 个交叉项（Cross_T1S1 ... Cross_T7S1）
  - 输出：trend_Cross, sig_Cross, amp_Cross
  - 保存到 Trend_Results / Amplitude_Results

关键说明
- 时间向量 Time_Vec 必须为 1×N 行向量，符合 gmt_harmonic_new 的输入要求。
- 该阶段仅对已有数据做回归计算，不产生新的物理量。

依赖
- gmt_harmonic_new（D:\work\Final_version\01_base）
- 02_calculate 产出的 MAT 文件
