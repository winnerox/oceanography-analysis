%% Step3_Calc_Generic_HighOrder.m
% =========================================================================
% 功能：通用高阶混合项计算器 (支持任意阶数)
% 引擎：gsw_rho_mixed_derivatives (迭代DP版)
% 
% 策略：
%   1. 使用循环自动生成并求和所有混合项 (Cross Terms)。
%   2. 允许 Avg 和 Std 使用不同的阶数 (Order Mismatch is OK)。
% 
% 配置：
%   Avg State: 2 -> 3 阶 (收敛快)
%   Std State: 2 -> 8 阶 (收敛慢，需对齐纯项)
% =========================================================================
clear; clc;

%% 1. 配置区域
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
SaveName  = fullfile(OutputDir, 'EN4_CrossTerms_Order8_Mismatch.mat');

Years = 2005:2024; 
rho0 = 1035.0; 
Std_T = 0.0; Std_S = 35.0;

% === 核心配置：阶数设置 ===
MaxOrder_Avg = 3;  % 平均态算到 3 阶
MaxOrder_Std = 8;  % 标准态算到 8 阶！

if isempty(which('gsw_rho_mixed_derivatives'))
    error('❌ 缺少计算引擎 gsw_rho_mixed_derivatives.m');
end

%% 2. 初始化网格
fprintf('>> [1/5] 初始化网格...\n');
% 自动找文件
FList = dir(fullfile(DataDir, 'EN.4.2.2.analyses.c13.2005', '*.nc'));
SampleFile = fullfile(FList(1).folder, FList(1).name);
Lon = ncread(SampleFile, 'lon');
Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= 2000);
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

%% 3. 准备参考态与预计算导数
fprintf('>> [2/5] 准备参考态与导数...\n');

% --- A. 平均态 (Avg) ---
fprintf('   计算 20 年平均场...\n');
Sum_T = zeros(Nx, Ny, Nz); Sum_S = zeros(Nx, Ny, Nz); Count = zeros(Nx, Ny, Nz, 'uint16');
for y = Years
    FPath = dir(fullfile(DataDir, sprintf('*%d*',y), '*.nc'));
    for k=1:length(FPath)
        try
            FF = fullfile(FPath(k).folder, FPath(k).name);
            T = double(ncread(FF,'temperature',[1 1 1 1],[Inf Inf Nz 1]));
            S = double(ncread(FF,'salinity',[1 1 1 1],[Inf Inf Nz 1]));
            if nanmean(T(:))>100, T=T-273.15; end
            Mask = ~isnan(T) & ~isnan(S);
            Sum_T(Mask)=Sum_T(Mask)+T(Mask); Sum_S(Mask)=Sum_S(Mask)+S(Mask); Count(Mask)=Count(Mask)+1;
        catch; end
    end
end
Mean_T = Sum_T./double(Count); Mean_T(Count==0)=NaN;
Mean_S = Sum_S./double(Count); Mean_S(Count==0)=NaN;
SA_Mean = gsw_SA_from_SP(Mean_S, P_3D, LON_3D, LAT_3D);
CT_Mean = gsw_CT_from_pt(SA_Mean, Mean_T);

% --- B. 标准态 (Std) ---
SA_Std_3D = Std_S * ones(Nx, Ny, Nz);
CT_Std_3D = Std_T * ones(Nx, Ny, Nz);

% --- C. 预计算所有需要的导数 ---
% 使用 Cell 数组存储导数矩阵，避免几百个变量名
% Deriv_Avg{n, k} 存储 ∂^n ρ / ∂S^k ∂T^(n-k)
fprintf('   预计算导数 (Avg: Max %d, Std: Max %d)...\n', MaxOrder_Avg, MaxOrder_Std);

% 预分配 Cell
Derivs_Avg = cell(MaxOrder_Avg, MaxOrder_Avg); 
Derivs_Std = cell(MaxOrder_Std, MaxOrder_Std);

% 计算 Avg 导数
for n = 2:MaxOrder_Avg
    for k = 1:n-1 % 混合项：k从1到n-1 (纯项是k=0和k=n，这里不存)
        fprintf('     Avg Deriv: Order %d (S^%d T^%d)\n', n, k, n-k);
        Derivs_Avg{n, k} = reshape(gsw_rho_mixed_derivatives.derivative(SA_Mean, CT_Mean, P_3D, k, n-k, 'S', 'T'), Nx, Ny, Nz);
    end
end

% 计算 Std 导数
for n = 2:MaxOrder_Std
    for k = 1:n-1
        fprintf('     Std Deriv: Order %d (S^%d T^%d)\n', n, k, n-k);
        Derivs_Std{n, k} = reshape(gsw_rho_mixed_derivatives.derivative(SA_Std_3D, CT_Std_3D, P_3D, k, n-k, 'S', 'T'), Nx, Ny, Nz);
    end
end

%% 4. 遍历计算
fprintf('>> [3/5] 开始遍历计算 (自动循环求和)...\n');
TotalFiles = length(Years)*12;
Cross_Avg_Total = zeros(Nx, Ny, TotalFiles, 'single');
Cross_Std_Total = zeros(Nx, Ny, TotalFiles, 'single');
time_vec = zeros(TotalFiles, 1);

% 预先计算阶乘倒数，减少循环内计算
InvFact = zeros(20,1); for i=1:20, InvFact(i)=1/factorial(i); end

idx = 0;
for y = Years
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        
        % 读取数据
        FPath = dir(fullfile(DataDir, sprintf('*%d*',y), sprintf('*%d%02d*.nc',y,m)));
        if isempty(FPath), continue; end
        FullPath = fullfile(FPath(1).folder, FPath(1).name);
        try
            T = double(ncread(FullPath,'temperature',[1 1 1 1],[Inf Inf Nz 1]));
            S = double(ncread(FullPath,'salinity',[1 1 1 1],[Inf Inf Nz 1]));
        catch; continue; end
        if nanmean(T(:))>100, T=T-273.15; end
        
        SA = gsw_SA_from_SP(S, P_3D, LON_3D, LAT_3D);
        CT = gsw_CT_from_pt(SA, T);
        
        % === Track A: Avg State (Loop up to MaxOrder_Avg) ===
        dS = SA - SA_Mean; dT = CT - CT_Mean;
        rho_sum_avg = zeros(Nx, Ny, Nz);
        
        for n = 2:MaxOrder_Avg
            for k = 1:n-1 % 混合项循环
                % Coeff = 1 / (k! * (n-k)!)
                C = InvFact(k) * InvFact(n-k);
                % Term = C * Deriv * dS^k * dT^(n-k)
                % 使用 bsxfun 确保幂运算安全，或者直接 .^ (如果版本够新)
                term = C * Derivs_Avg{n,k} .* (dS.^k) .* (dT.^(n-k));
                rho_sum_avg = rho_sum_avg + term;
            end
        end
        Cross_Avg_Total(:,:,idx) = -(nansum(rho_sum_avg .* dz_perm, 3) / rho0) * 1000;
        
        % === Track B: Std State (Loop up to MaxOrder_Std) ===
        dS = SA - Std_S; dT = CT - Std_T;
        rho_sum_std = zeros(Nx, Ny, Nz);
        
        for n = 2:MaxOrder_Std
            for k = 1:n-1
                C = InvFact(k) * InvFact(n-k);
                term = C * Derivs_Std{n,k} .* (dS.^k) .* (dT.^(n-k));
                rho_sum_std = rho_sum_std + term;
            end
        end
        Cross_Std_Total(:,:,idx) = -(nansum(rho_sum_std .* dz_perm, 3) / rho0) * 1000;
        
        if mod(idx, 24)==0, fprintf('   Progress: %.1f%%\n', idx/TotalFiles*100); end
    end
end

%% 5. 保存
fprintf('>> [4/5] 保存结果...\n');
% 保持变量名兼容性，以便后续绘图脚本直接用
HO_Cross_Avg = Cross_Avg_Total;
HO_Cross_Std_TSLA = Cross_Std_Total; % 注意这里是 Total

save(SaveName, 'HO_Cross_Avg', 'HO_Cross_Std_TSLA', 'Lon', 'Lat', 'time_vec', '-v7.3');
fprintf('>> [5/5] ✅ 完成！\n');
fprintf('   Avg 阶数: %d\n', MaxOrder_Avg);
fprintf('   Std 阶数: %d\n', MaxOrder_Std);
fprintf('   现在请修改绘图脚本，加载这个新的 .mat 文件。\n');