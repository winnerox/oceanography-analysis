%% Calc_S_Terms_IAP_Avg.m
% =========================================================================
% 功能：计算 IAP【平均态】下的纯盐度项 (Pure S Terms, 1-8阶)
% 输出：D:\work\MAT_Data\IAP_S_Terms_1to8_Average.mat
% 极致优化：
%   1. 掩码逻辑：只对有效深度求和，保留浅海数据
%   2. 强制单精度：所有变量显式single类型
% =========================================================================
clear; clc;
addpath('D:\work');

DataDir = 'D:\work\IAP_05_24';
OutputDir = 'D:\work\MAT_Data';
MeanFile = 'D:\work\IAP_TSLA_Terms\IAP_Mean_State.mat';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'IAP_S_Terms_1to8_Average.mat');

Years = 2005:2024; rho0 = single(1035.0); MaxDepth = 2000; MaxOrder = 8;
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

if ~exist(MeanFile, 'file'), error('请先运行 Calc_Mean_IAP.m'); end
load(MeanFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');
Mean_T = single(Mean_T); Mean_S = single(Mean_S);
Lon = single(Lon); Lat = single(Lat);
CalcDepth = single(Depth(DepthIdx));

[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

fprintf('Calculating Derivatives...\n');
SA_Ref = single(gsw_SA_from_SP(double(Mean_S), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref = single(gsw_CT_from_pt(SA_Ref, double(Mean_T)));
Derivs_S = cell(MaxOrder, 1);
for n = 1:MaxOrder
    D_val = Engine.calculate_mixed(SA_Ref, CT_Ref, P_3D, 0, n);
    Derivs_S{n} = single(D_val);
end

InvFact = zeros(21,1,'single'); for i=0:20, InvFact(i+1)=single(1/factorial(i)); end

TotalSteps = length(Years) * 12;
SSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_vec = zeros(TotalSteps, 1);
idx = 0;

for y = Years
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        Pat = sprintf('*year_%d_month_%02d.nc', y, m);
        
        try
            sS = dir(fullfile(DataDir, 'SALT', Pat));
            if isempty(sS), continue; end
            FS = fullfile(sS(1).folder, sS(1).name);
            
            S_raw = single(ncread(FS, 'salinity', [1 1 1], [Inf Inf Nz]));
            
            MaskNaN = isnan(S_raw) | (S_raw < 0) | (S_raw >= 100);
            S_raw(MaskNaN) = 0;
            
            SA_Inst = single(gsw_SA_from_SP(double(S_raw), double(P_3D), double(LON_3D), double(LAT_3D)));
            dS = SA_Inst - SA_Ref;
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
