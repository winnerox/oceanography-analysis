import h5py
import os
import sys

files = [
    'D:/work/IAP_05_24/TEMP/IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc',
    'D:/work/Ishii_05_24/Salinity/sal.2005.nc'
]

for f in files:
    print(f"--- {f} ---")
    try:
        with h5py.File(f, 'r') as hf:
            print("HDF5 Keys:", list(hf.keys()))
            for k in hf.keys():
                try: print(f"  {k}: {hf[k].shape}")
                except: pass
    except Exception as e:
        print("HDF5 failed:", e)
        try:
            from scipy.io import netcdf
            with netcdf.netcdf_file(f, 'r') as nf:
                print("NetCDF3 Keys:", list(nf.variables.keys()))
                for k in nf.variables.keys():
                    try: print(f"  {k}: {nf.variables[k].shape}")
                    except: pass
        except Exception as e2:
            print("NetCDF3 failed:", e2)
