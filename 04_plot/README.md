# 04_plot

Purpose
- Generate spatial maps, boxplots, ensemble statistics, and diagnostics from trend/amplitude products.

Inputs
- D:\work\Task_Convergence\Trend_Results\*
- D:\work\Task_Convergence\Amplitude_Results\*

Outputs
- Figures saved under D:\work\Figures (some scripts also display interactively).

Key scripts
- D:\work\Final_version\04_plot\Ensemble_Analysis_AllInOne.m
  - Multi-dataset ensemble mean and confidence intervals; optional spatial interpolation to a common grid.
- D:\work\Final_version\04_plot\Plot_Trend_Combined.m
  - Spatial trend maps and boxplots for SSLA/TSLA/HSLA, with significance masking.
- D:\work\Final_version\04_plot\Plot_Amplitude_Combined.m
  - Spatial amplitude maps and boxplots for SSLA/TSLA/HSLA.
- D:\work\Final_version\04_plot\Plot_CrossDetail_Trend.m
  - Trend maps for cross terms (TnSm).
- D:\work\Final_version\04_plot\Plot_Cumulative_Sum.m
  - Cumulative sum diagnostics.
- D:\work\Final_version\04_plot\Plot_Exact_SSLA_Difference.m
  - Compare exact SSLA between datasets/states.
- D:\work\Final_version\04_plot\Plot_State_Difference.m
  - Compare Average vs StdRef states.
- D:\work\Final_version\04_plot\Plot_Steric_Budget.m
  - Steric budget decomposition plots.
- D:\work\Final_version\04_plot\Plot_Trend_Combined.m
  - Combined trend visualization with adaptive scaling and grouping.

Notes
- Most scripts expect Lon/Lat fields named consistently (Lon, Lat) from prior steps.
- Output filenames are built from dataset/state/variable names.
