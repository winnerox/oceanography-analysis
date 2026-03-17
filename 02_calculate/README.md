# 02_calculate

Purpose
- Compute steric sea level components for three datasets with two reference states.
- Produces time series grids for SSLA (exact), TSLA/HSLA Taylor terms (1-8), and Cross terms.

Inputs (local data)
- EN4: D:\work\EN4_analyses_c13_last20years
- IAP: D:\work\IAP_05_24 (TEMP, SALT)
- Ishii: D:\work\Ishii_05_24 (Temperature, Salinity)

Outputs (MAT files)
- D:\work\EN4_mat_data\*
- D:\work\IAP_mat_data\*
- D:\work\Ishii_mat_data\*

Key scripts
- D:\work\Final_version\02_calculate\EN4_4in1.m
  - Builds 2005-2024 mean state, computes TEOS-10 derivatives, and outputs SSLA/TSLA/HSLA/Cross for Average and StdRef states.
- D:\work\Final_version\02_calculate\IAP_4in1.m
  - Same workflow for IAP with salinity longitude interpolation and depth-bounds reconstruction.
- D:\work\Final_version\02_calculate\Ishii_4in1.m
  - Same workflow for Ishii yearly files (monthly slices).

Notes
- Uses TEOS10_General_Engine and GSW functions.
- Max depth default is 2000 m; max Taylor order default is 8.
