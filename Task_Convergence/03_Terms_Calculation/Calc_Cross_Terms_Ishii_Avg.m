%% Calc_Cross_Terms_Ishii_Avg.m
% =========================================================================
% 功能：计算 Ishii【平均态】下的混合项 (Cross Terms, 2-8阶)
% 输出：D:\work\MAT_Data\Ishii_Cross_Terms_1to8_Average.mat
% 极致优化：
%   1. 避免4D内存杀手：不预计算S_Pow/T_Pow，现场计算幂次
%   2. 强制单精度：所有变量显式single类型
%   3. 一阶冗余：从n=2开始
% =========================================================================
clear; clc;
addpath('D:\work');

DataDir = 'D:\work\Ishii_05_24';
OutputDir = 'D:\work\MAT_Data';
MeanStateFile = 'D:\work\Ishii_TSLA_Terms\Ishii_Mean_State.mat';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'Ishii_Cross_Terms_1to8_Average.mat');

Years = 2005:2024; rho0 = single(1035.0); MaxDepth = 2000; MaxOrder = 8;

fprintf('[1/6] 加载平均态数据...\n');
load(MeanStateFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');
Mean_T = single(Mean_T); Mean_S = single(Mean_S);
Lon = single(Lon); Lat = single(Lat);
CalcDepth = single(Depth(DepthIdx));

fprintf('[2/6] 初始化网格...\n');
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

fprintf('[3/6] 初始化 TEOS-10 引擎...\n');
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

SA_Ref = single(gsw_SA_from_SP(double(Mean_S), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref = single(gsw_CT_from_pt(SA_Ref, double(Mean_T)));

Derivs = cell(MaxOrder, MaxOrder+1);
for n = 1:MaxOrder
    for k = 0:n
        n_T = n - k; n_S = k;
        D_val = Engine.calculate_mixed(SA_Ref, CT_Ref, P_3D, n_T, n_S);
        Derivs{n, k+1} = single(D_val);
    end
end

InvFact = zeros(21,1,'single'); for i=0:20, InvFact(i+1)=single(1/factorial(i)); end

TotalSteps = length(Years) * 12;
Cross_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_vec = zeros(TotalSteps, 1);
idx = 0;

for y = Years
    FT = fullfile(DataDir, 'Temperature', sprintf('temp.%d.nc', y));
    FS = fullfile(DataDir, 'Salinity', sprintf('sal.%d.nc', y));
    
    if ~exist(FT, 'file') || ~exist(FS, 'file')
        for m = 1:12
            idx = idx + 1;
            time_vec(idx) = y + (m-0.5)/12;
            Cross_AllOrders(:,:,idx,:) = NaN;
        end
        continue;
    end
    
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        
        try
            T_raw = single(ncread(FT, 'temp', [1 1 1 m], [Inf Inf Nz 1]));
            S_raw = single(ncread(FS, 'sal', [1 1 1 m], [Inf Inf Nz 1]));
            
            MaskNaN = isnan(T_raw) | isnan(S_raw) | (T_raw <= -100) | (T_raw >= 100) | (S_raw < 0) | (S_raw >= 100);
            T_raw(MaskNaN) = 0; S_raw(MaskNaN) = 0;
            
            SA_Inst = single(gsw_SA_from_SP(double(S_raw), double(P_3D), double(LON_3D), double(LAT_3D)));
            CT_Inst = single(gsw_CT_from_pt(SA_Inst, double(T_raw)));
            
            dS = SA_Inst - SA_Ref;
            dT = CT_Inst - CT_Ref;
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
