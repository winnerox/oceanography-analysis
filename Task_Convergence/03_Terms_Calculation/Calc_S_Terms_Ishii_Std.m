%% Calc_S_Terms_Ishii_Std.m
% =========================================================================
% 功能：计算 Ishii【标准态】下的纯盐度项 (Pure S Terms, 1-8阶)
% 输出：D:\work\MAT_Data\Ishii_S_Terms_1to8_StdRef.mat
% 极致优化：
%   1. 掩码逻辑：只对有效深度求和，保留浅海数据
%   2. 强制单精度：所有变量显式single类型
% =========================================================================
clear; clc;
addpath('D:\work');

DataDir = 'D:\work\Ishii_05_24';
OutputDir = 'D:\work\MAT_Data';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'Ishii_S_Terms_1to8_StdRef.mat');

Years = 2005:2024; rho0 = single(1035.0); MaxDepth = 2000; MaxOrder = 8;
Std_S_Val = single(35.0); Std_T_Val = single(0.0);
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

List = dir(fullfile(DataDir, 'Temperature', '*.nc'));
SampleFile = fullfile(List(1).folder, List(1).name);

vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'lon')), vLon='lon'; elseif any(strcmp(vnames, 'long')), vLon='long'; else, vLon=vnames{1}; end
if any(strcmp(vnames, 'lat')), vLat='lat'; elseif any(strcmp(vnames, 'latitude')), vLat='latitude'; else, vLat=vnames{2}; end
if any(strcmp(vnames, 'level')), vDep='level';
elseif any(strcmp(vnames, 'depth')), vDep='depth'; 
elseif any(strcmp(vnames, 'depth_std')), vDep='depth_std'; 
else, vDep='level'; end

Lon = single(ncread(SampleFile, vLon)); Lat = single(ncread(SampleFile, vLat));
Depth = single(ncread(SampleFile, vDep)); 

if any(isnan(Lon)), Lon(isnan(Lon)) = 0; end
if any(isnan(Lat)), Lat(isnan(Lat)) = 0; end
Lon = single(mod(Lon, 360)); 
Lat(Lat < -90) = single(-90); Lat(Lat > 90) = single(90);

DepthIdx = find(Depth <= MaxDepth); 
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

fprintf('Calculating Derivatives...\n');

% 使用3D参考场，而非1D捷径！(与空间坐标保持一致)
SA_Ref_3D = single(gsw_SA_from_SP(double(Std_S_Val * ones(Nx, Ny, Nz, 'single')), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref_3D = single(gsw_CT_from_pt(SA_Ref_3D, double(Std_T_Val * ones(Nx, Ny, Nz, 'single'))));

Derivs_S = cell(MaxOrder, 1);
for n = 1:MaxOrder
    D_val = Engine.calculate_mixed(double(SA_Ref_3D), double(CT_Ref_3D), double(P_3D), 0, n);
    Derivs_S{n} = single(D_val);
end

InvFact = zeros(21,1,'single'); for i=0:20, InvFact(i+1)=single(1/factorial(i)); end

TotalSteps = length(Years) * 12;
SSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_vec = zeros(TotalSteps, 1);
idx = 0;

for y = Years
    FS = fullfile(DataDir, 'Salinity', sprintf('sal.%d.nc', y));
    if ~exist(FS, 'file')
        for m = 1:12
            idx = idx + 1;
            time_vec(idx) = y + (m-0.5)/12;
            SSLA_AllOrders(:,:,idx,:) = NaN;
        end
        continue;
    end
    
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        
        try
            S_raw = single(ncread(FS, 'sal', [1 1 1 m], [Inf Inf Nz 1]));
            S_raw(S_raw > 100) = NaN;
            
            MaskNaN = isnan(S_raw) | (S_raw < 0) | (S_raw >= 100);
            S_raw(MaskNaN) = 0;
            
            SA_Inst = single(gsw_SA_from_SP(double(S_raw), double(P_3D), double(LON_3D), double(LAT_3D)));
            
            dS = SA_Inst - SA_Ref_3D;
            dS(MaskNaN) = 0;
            
            for n = 1:MaxOrder
                C = InvFact(n+1);
                Deriv = Derivs_S{n};
                
                rho_s = C .* Deriv .* (dS .^ n);
                rho_s(MaskNaN) = 0;
                
                val = -single(sum(rho_s .* dz_perm, 3)) / rho0 * single(1000);
                
                valid_depths = single(sum(~MaskNaN, 3));
                val(valid_depths < 1) = NaN;
                
                SSLA_AllOrders(:,:,idx,n) = val;
            end
        catch
            SSLA_AllOrders(:,:,idx,:) = NaN;
        end
    end
    fprintf('Year %d Done.\n', y);
end

lon = Lon; lat = Lat;
save(SaveName, 'SSLA_AllOrders', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('Saved %s\n', SaveName);
