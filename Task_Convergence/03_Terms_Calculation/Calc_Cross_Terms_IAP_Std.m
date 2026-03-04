%% Calc_Cross_Terms_IAP_Std.m
% =========================================================================
% 功能：计算 IAP【标准态】下的混合项 (Cross Terms, 2-8阶)
% 输出：D:\work\MAT_Data\IAP_Cross_Terms_1to8_StdRef.mat
% 极致优化：
%   1. 避免4D内存杀手：不预计算S_Pow/T_Pow，现场计算幂次
%   2. 强制单精度：所有变量显式single类型
%   3. 一阶冗余：从n=2开始
% =========================================================================
clear; clc;
addpath('D:\work');

DataDir = 'D:\work\IAP_05_24';
OutputDir = 'D:\work\MAT_Data';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'IAP_Cross_Terms_1to8_StdRef.mat');

Years = 2005:2024; rho0 = single(1035.0); MaxDepth = 2000; MaxOrder = 8;
Std_S_Val = single(35.0); Std_T_Val = single(0.0);
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

List = dir(fullfile(DataDir, 'TEMP', '*.nc'));
SampleFile = fullfile(List(1).folder, List(1).name);

vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth_std')), vDepth='depth_std';
elseif any(strcmp(vnames, 'depth')), vDepth='depth'; 
else, vDepth='depth_std'; end

Lon = single(ncread(SampleFile, 'lon')); Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, vDepth)); 
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
Std_S = ones(Nx, Ny, Nz, 'single') * Std_S_Val;
Std_T = ones(Nx, Ny, Nz, 'single') * Std_T_Val;

% 使用3D参考场，而非1D捷径！(与空间坐标保持一致)
SA_Ref_3D = single(gsw_SA_from_SP(double(Std_S), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref_3D = single(gsw_CT_from_pt(SA_Ref_3D, double(Std_T)));

Derivs = cell(MaxOrder, MaxOrder+1);
for n = 1:MaxOrder
    for k = 0:n
        n_T = n - k; n_S = k;
        D_val = Engine.calculate_mixed(double(SA_Ref_3D), double(CT_Ref_3D), double(P_3D), n_T, n_S);
        Derivs{n, k+1} = single(D_val);
    end
end

InvFact = zeros(21,1,'single'); for i=0:20, InvFact(i+1)=single(1/factorial(i)); end

TotalSteps = length(Years) * 12;
Cross_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_vec = zeros(TotalSteps, 1);
idx = 0;

for y = Years
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        Pat = sprintf('*year_%d_month_%02d.nc', y, m);
        
        try
            sT = dir(fullfile(DataDir, 'TEMP', Pat));
            sS = dir(fullfile(DataDir, 'SALT', Pat));
            if isempty(sT) || isempty(sS), continue; end
            
            FT = fullfile(sT(1).folder, sT(1).name);
            FS = fullfile(sS(1).folder, sS(1).name);
            
            T_raw = single(ncread(FT, 'temp', [1 1 1], [Inf Inf Nz]));
            S_raw = single(ncread(FS, 'salinity', [1 1 1], [Inf Inf Nz]));
            
            MaskNaN = isnan(T_raw) | isnan(S_raw) | (T_raw <= -100) | (T_raw >= 100) | (S_raw < 0) | (S_raw >= 100);
            T_raw(MaskNaN) = 0; S_raw(MaskNaN) = 0;
            
            SA_Inst = single(gsw_SA_from_SP(double(S_raw), double(P_3D), double(LON_3D), double(LAT_3D)));
            CT_Inst = single(gsw_CT_from_pt(SA_Inst, double(T_raw)));
            
            dS = SA_Inst - SA_Ref_3D;
            dT = CT_Inst - CT_Ref_3D;
            dS(MaskNaN) = 0; dT(MaskNaN) = 0;
            
            for n = 2:MaxOrder
                rho_cross = zeros(Nx, Ny, Nz, 'single');
                for k = 1:n-1
                    C = InvFact(k+1) * InvFact(n-k+1);
                    Deriv = Derivs{n, k+1};
                    
                    term = C .* Deriv .* (dS .^ k) .* (dT .^ (n-k));
                    term(MaskNaN) = 0;
                    rho_cross = rho_cross + term;
                end
                val = -single(sum(rho_cross .* dz_perm, 3)) / rho0 * single(1000);
                
                valid_depths = single(sum(~MaskNaN, 3));
                val(valid_depths < 1) = NaN;
                
                Cross_AllOrders(:,:,idx,n) = val;
            end
        catch
            Cross_AllOrders(:,:,idx,:) = NaN;
        end
    end
    fprintf('Year %d Done.\n', y);
end

lon = Lon; lat = Lat;
save(SaveName, 'Cross_AllOrders', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('Saved %s\n', SaveName);
