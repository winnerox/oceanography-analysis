# Copilot Instructions - Ocean Temperature & Salinity Analysis

## Project Overview
This project analyzes **ocean steric sea-level change** using EN4 oceanographic data. It computes high-order derivatives of density with respect to temperature and salinity using TEOS-10 thermodynamic equations, then performs trend analysis and nonlinearity diagnosis on the resulting **Thermosteric Sea-Level Anomaly (TSLA)** terms across 1st-8th order.

### Key Domains
- **MATLAB**: Main analysis pipeline, visualization, statistical computations
- **Python**: Data downloading, TEOS-10 physics calculations (dual implementations)
- **Data Format**: NetCDF (EN4 oceanic profiles), MAT files (cached results)
- **Physics**: TEOS-10 equation of state, Leibniz rule for derivative composition, thermodynamic properties (density, specific volume)

---

## Architecture & Data Flow

### 1. **Core Component: TEOS-10 Derivative Engine**
Both MATLAB (`TEOS10_General_Engine.m`) and Python (`teos10_engine.py`) implement identical **Leibniz-based recursive differentiation**:

```
Input: (SA, CT, p) → Output: High-order partial derivatives
drho/dT^n, drho/dS^n at orders 1-8
```

**Key Implementation Detail**: GSW (Gibbs SeaWater) coefficients stored in `gsw_coeffs.txt` drive polynomial expansion. Classes precompute binomial coefficients to avoid recursion overhead.

**MATLAB**: Vectorized matrix operations in `calc_leibniz_vectorized()` (see line ~120 in `TEOS10_General_Engine.m`)  
**Python**: Numpy array-based computation with optional GPU acceleration via CuPy

---

### 2. **Data Processing Pipeline**
```
EN4 Raw Data (.nc) 
   ↓ [download_IAP_lshii.py]
   ↓ [Retry logic, parallel wget with ThreadPoolExecutor]
   ↓ Cache: Ishii_v7.3.1_Data_05_24/
   
EN4 Data Processing (2005-2024)
   ↓ [Step1_Calc_TSLA_HighOrder.m OR step1_calc_python.py]
   ↓ Load each monthly file, extract (SA, CT, p) at grid points
   ↓ Apply TEOS-10 engine to compute ρ derivatives
   ↓ Convert to TSLA via: ΔSL = -∫Δρ/ρ₀ dz (integrated over depth 0-2000m)
   ↓ Output: EN4_TSLA_Terms1to8_Average.mat
   
TSLA Analysis & Visualization
   ↓ [All_in_One_Analysis_EN4.m]
   ↓ Compute trend (mm/yr), amplitude (std), statistical significance
   ↓ Generate 2D maps (space-time) + box plots
   ↓ Output figures to output/
```

---

### 3. **Reference States & Standardization**
- **Average Anomaly**: Reference is 20-year mean state (rolling climatology)
- **Standard Reference**: Fixed standard seawater (S=35, T=0°C)
  - Files: `*_StdRef.mat` variant scripts use `Step1_Calc_TSLA_StdRef_RealSalt.m`
  - Purpose: Compare nonlinearity impact vs. fixed reference

**Critical Convention**: All depths in metadata use **positive values** (0-2000 m); code converts to negative when calling GSW functions via `gsw.p_from_z(-depth, lat)` (see `step1_calc_python.py` line 39).

---

## File Organization Patterns

### Numbered Step Scripts
Scripts follow naming: `Step{NN}_{Task}.m` (primary) or corresponding Python variant:
- **Step1**: Calculate TSLA derivatives and store base matrices
- **Step18-24**: Spatial analysis, derivatives, mechanism diagnosis
- **Step31**: Compute nonlinearity index (ratio of higher-order corrections)

**Design Pattern**: Each step is idempotent—can rerun independently given prerequisite MAT files.

### Output Caching Strategy
Large computations cached in `.mat` files to avoid recomputation:
- `EN4_TSLA_Terms1to8_Average.mat` → Contains `TSLA_AllOrders [nlon×nlat×ntime×8]`
- `Trend_Result_Final.mat` → Contains `Map_Trend, Map_Sig, Map_Amp [nlat×nlon×8]`
- `EN4_TSLA_Nonlinearity_Index.mat` → Stores computed nonlinearity metrics

**Load Pattern** (see `All_in_One_Analysis_EN4.m` line 15-30):
```matlab
load('EN4_TSLA_Terms/EN4_TSLA_Terms1to8_Average.mat')
% Expect: TSLA_AllOrders, time_axis, lat, lon
```

---

## Critical Coding Conventions

### 1. **Dimension Ordering**
**MATLAB** (column-major):
- Data indexing: `(longitude, latitude, time, order)` → `[360, 173, 240, 8]`
- Output maps: **transpose to** `(latitude, longitude, order)` → `[173, 360, 8]` for storage
- **Common Bug**: Forgetting this transpose when assigning to result matrices

Example Fix (from `All_in_One_Analysis_EN4.m` line 58):
```matlab
temp_trend = nan(nlat, nlon);  % Pre-allocate as [nlat, nlon] NOT [nlon, nlat]
for i = 1:nlon
    for j = 1:nlat
        ... = temp_trend(j, i);  % Access as (lat_idx, lon_idx)
    end
end
```

### 2. **Longitude Wrapping**
EN4 data uses 0-360°. Code normalizes to -180-180°:
```matlab
if max(lon) > 180
    lon(lon > 180) = lon(lon > 180) - 360;
end
[lon, sort_idx] = sort(lon);  % MUST re-sort after wrapping
TSLA_AllOrders = TSLA_AllOrders(sort_idx, :, :, :);  % Permute data accordingly
```

### 3. **Unit Conventions**
- **TSLA**: millimeters (mm), computed as `ΔSL × 1000`
- **Trend**: mm/year (divide time in decades by result)
- **Depth**: positive values in metadata (0-2000), converted to negative for GSW calls
- **Pressure**: dbar, computed via `gsw.p_from_z(-depth, latitude)`
- **Practical Salinity (SA)**: PSU; **Absolute Salinity (SA)** ≠ practical—use GSW conversion if switching

---

## Physics & Mathematical Notes

### TEOS-10 Derivatives as Taylor Series
The engine implements density as a polynomial in scaled variables (Xs, Xt, Xp):
```
ρ(SA, CT, p) = Σ aijk·(Xs)^i·(Xt)^j·(Xp)^k
```
where Xs = SA + 24, Xt = CT, Xp = p, with standard scalings `[0.0248826, 0.025, 1e-4]`.

**Leibniz Rule for Derivatives**:
- `dⁿρ/dT^n` computed via Leibniz composition: `dⁿρ = -1/ρ₀ · Σ C(n,k)·dᵏρ·dⁿ⁻ᵏv`
- Avoids explicit high-order differentiation of polynomial; instead chains first-order terms

### Nonlinearity Index
Computed as ratio of cumulative high-order corrections to linear term:
```
Nonlinearity_Index = |Σ(d²-d⁸ terms)| / |d¹ term|
```
See `Step31_Calc_Nonlinearity_Index.m` for implementation.

---

## Data Dependencies & External Tools

### Required Data Sources
1. **EN4 Oceanographic Profiles**: https://climate.mri-jma.go.jp/pub/ocean/ts/v7.3.1/
   - Downloaded via `download_IAP_lshii.py` with retry logic & parallel workers
   - Expected structure: `EN4_analyses_c13_last20years/{YEAR}/{MONTH}_files.nc`

2. **GSW Matlab Toolbox**: `/gsw_matlab_v3_06_16/` folder must be in path
   - Alternative: Python's `gsw` package (via pip)

3. **CDT (Climate Data Toolbox)**: `/CDT-master/` for visualization helpers
   - Used in mapping scripts (`Step20_Map_Trend_Final.m`) for globe-level plots

### Configuration Hardcoding
User paths hardcoded in scripts; common variables to update:
```matlab
% Top of most Step*.m files:
DataDir = 'D:\work\EN4_analyses_c13_last20years';  % ← Update to local EN4 path
OutputDir = 'D:\work\EN4_TSLA_Terms';              % ← Output cache folder
Years = 2005:2024;
MaxDepth = 2000;  % Depth cutoff (m)
rho0 = 1028.0;    % Reference density (kg/m³)
```

---

## Common Workflows & Debugging

### Full Recomputation (Rare)
1. Ensure EN4 data exists: `D:\work\EN4_analyses_c13_last20years/`
2. Run `Step1_Calc_TSLA_HighOrder.m` → generates `EN4_TSLA_Terms1to8_Average.mat`
3. Run trend analysis: `All_in_One_Analysis_EN4.m`
4. Visualize: `Step20_Map_Trend_Final.m`, `Step21_Map_Amp_Final.m`

**Expected Runtime**: ~1-2 hrs (MATLAB) or ~30 min (Python + GPU)

### Quick Diagnosis
- **Dimension mismatch errors**: Check file load expectations against actual MAT structure (use `whos` in MATLAB)
- **NaN propagation**: Enable `Check for NaN/Inf` in MATLAB settings; validate input EN4 files contain valid (SA, CT, p)
- **Trend anomalies**: Visually inspect `test_TEOS10_convergence.fig` to confirm engine convergence at order 8

### Python↔MATLAB Compatibility
Both implement identical TEOS-10, but outputs differ in rounding at ~12th decimal due to numerical precision. Validate cross-checks via unit tests (see `test_TEOS10_General_Engine.m`).

---

## Key Files Reference
| File | Purpose |
|------|---------|
| `TEOS10_General_Engine.m` | Core MATLAB physics engine |
| `teos10_engine.py` | Core Python physics engine (CuPy support) |
| `Step1_Calc_TSLA_HighOrder.m` | Primary TSLA computation pipeline |
| `All_in_One_Analysis_EN4.m` | Trend/amplitude analysis & visualization |
| `download_IAP_lshii.py` | Data acquisition with retry logic |
| `gsw_matlab_v3_06_16/` | Gibbs SeaWater reference implementation |
| `EN4_TSLA_Terms/` | Output cache (preprocessed matrices) |
