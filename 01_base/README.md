# 01_base

Purpose
- Shared core math/physics utilities used by the pipeline.

Key files
- D:\work\Final_version\01_base\gmt_harmonic_new.m
  - Harmonic regression on a time series (trend + annual + semi-annual).
  - Returns amplitudes/phases, trend statistics (t, p, significance), residuals, and optional interpolated series.
  - Supports 2D/3D grids (lat, lon, time) with optional per-gridpoint std.
- D:\work\Final_version\01_base\TEOS10_General_Engine.m
  - TEOS-10 high-order derivative engine for seawater properties.
  - Vectorized computation of pure and mixed derivatives up to a chosen order.

Typical use
- Called by dataset pipelines in D:\work\Final_version\02_calculate to build steric components.
- Called by trend/amplitude builders in D:\work\Final_version\03_calculate_trend_amplitude.
