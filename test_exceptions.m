clear; clc;
addpath('D:\work');
try
    disp('--- Testing IAP Cross Terms ---');
    y = 2005; m = 1;
    DataDir = 'D:\work\IAP_05_24';
    Pat = sprintf('*year_%d_month_%02d.nc', y, m);
    sT = dir(fullfile(DataDir, 'TEMP', Pat));
    sS = dir(fullfile(DataDir, 'SALT', Pat));
    FT = fullfile(sT(1).folder, sT(1).name);
    FS = fullfile(sS(1).folder, sS(1).name);
    
    Nz = 41;
    T_raw = single(permute(ncread(FT, 'temp', [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
    S_raw = single(permute(ncread(FS, 'salinity', [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
    disp('IAP Read Success!');
catch ME
    disp('IAP Error:');
    disp(ME.message);
end

try
    disp('--- Testing Ishii S Terms ---');
    y = 2005; m = 1;
    DataDir = 'D:\work\Ishii_05_24';
    FS = fullfile(DataDir, 'Salinity', sprintf('sal.%d.nc', y));
    Nz = 24;
    S_raw = single(ncread(FS, 'sal', [1 1 1 m], [Inf Inf Nz 1]));
    disp('Ishii S Terms Read Success!');
catch ME
    disp('Ishii Error:');
    disp(ME.message);
end

% Check Calc_S_Terms_Ishii_Avg loop explicitly
try
    disp('--- Testing Ishii Derivs ---');
    MeanFile = 'D:\work\Ishii_TSLA_Terms\Ishii_Mean_State.mat';
    load(MeanFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');
    Mean_T = single(Mean_T); Mean_S = single(Mean_S);
    Lon = single(Lon); Lat = single(Lat);
    if any(isnan(Lon)), Lon(isnan(Lon)) = 0; end
    if any(isnan(Lat)), Lat(isnan(Lat)) = 0; end
    CalcDepth = single(Depth(DepthIdx));
    [LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
    P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
    SA_Ref = single(gsw_SA_from_SP(double(Mean_S), double(P_3D), double(LON_3D), double(LAT_3D)));
    CT_Ref = single(gsw_CT_from_pt(SA_Ref, double(Mean_T)));
    
    MaxOrder = 8;
    CacheFile = fullfile('D:\work\MAT_Data', 'TEOS10_Engine_Cache.mat');
    if exist(CacheFile, 'file')
        Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);
    else
        Engine = TEOS10_HighOrder_Engine(MaxOrder, '', false);
    end
    
    Derivs_S = cell(MaxOrder, 1);
    for n = 1:MaxOrder
        D_val = Engine.calculate_mixed(SA_Ref, CT_Ref, P_3D, 0, n);
        Derivs_S{n} = single(D_val);
    end
    disp('Ishii Derivs Success!');
catch ME
    disp('Ishii Derivs Error:');
    disp(ME.message);
end
exit;
