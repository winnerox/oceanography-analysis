# 03_calculate_trend_amplitude

Purpose
- Build trend and seasonal amplitude products from the steric time series.
- Uses harmonic regression (gmt_harmonic_new) for SSLA/TSLA/HSLA and Cross terms.

Inputs
- D:\work\EN4_mat_data\*
- D:\work\IAP_mat_data\*
- D:\work\Ishii_mat_data\*

Outputs
- D:\work\Task_Convergence\Trend_Results\*
- D:\work\Task_Convergence\Amplitude_Results\*

Key scripts
- D:\work\Final_version\03_calculate_trend_amplitude\Build_Trend_Amplitude_Data.m
  - Computes trend/significance and annual amplitude for SSLA exact + TSLA/HSLA orders 1-8.
- D:\work\Final_version\03_calculate_trend_amplitude\Build_Cross_Trend_Amplitude.m
  - Computes trend/significance and annual amplitude for 28 cross terms (T1S1 ... T7S1).

Notes
- Time vector is expected as row vector (1 x N) for gmt_harmonic_new.
- Outputs are structured for plotting scripts in D:\work\Final_version\04_plot.
