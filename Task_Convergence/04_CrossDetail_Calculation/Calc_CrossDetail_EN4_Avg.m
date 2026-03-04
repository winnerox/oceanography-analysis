%% Calc_CrossDetail_Avg.m
% =========================================================================
% 功能：计算【平均态】下的混合偏导数项（逐项分解）
% 阶数：2-5阶 (共 10 项)
% 引擎：TEOS10_HighOrder_Engine
% 
% 输出：EN4_CrossDetail_Avg.mat
%   包含变量如: Cross_T1S1, Cross_T2S1, Cross_T1S2, ...
% reference: Mean State (S_mean, T_mean)
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
MeanStateFile = fullfile(OutputDir, 'EN4_Mean_State.mat');
SaveName = fullfile(OutputDir, 'EN4_CrossDetail_Avg.mat');
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; 
rho0 = 1035.0;       
MaxDepth = 2000;     
MaxOrder = 8;

% 初始化引擎
CacheFile = fullfile(OutputDir, 'TEOS10_Engine_Cache.mat');
fprintf('>> [1/6] 初始化 TEOS10 高阶引擎 (MaxOrder=%d)...\n', MaxOrder);
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false); % Simplify=false

%% 2. 加载平均态
fprintf('>> [2/6] 加载平均态数据...\n');
if ~exist(MeanStateFile, 'file'), error('❌ 找不到平均态文件'); end
load(MeanStateFile, 'Mean_T', 'Mean_S'); 
Mean_T = double(Mean_T); Mean_S = double(Mean_S);

%% 3. 初始化网格
fprintf('>> [3/6] 初始化网格...\n');
FList = dir(fullfile(DataDir, 'EN.4.2.2.analyses.c13.2005', '*.nc'));
if isempty(FList), error('找不到数据文件'); end
SampleFile = fullfile(FList(1).folder, FList(1).name);

Lon = ncread(SampleFile, 'lon'); Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= MaxDepth); CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

% 垂直权重
dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

%% 4. 预计算平均态导数系数
fprintf('>> [4/6] 计算平均态导数系数 (1D Profile Approximation)...\n');
% 为节省内存，采用 Profile 预计算 (Mean State Spatial Variation kept for dS/dT, but Derivs approximated)
% 注意：平均态在空间上变化显著。
% 对于平均态分析，通常在当地平均态处计算导数。
% 之前的脚本 `Calc_Cross_and_S_Terms_Avg.m` 计算了全 3D 导数。
% 由于 MaxOrder=5 不算太高，且只计算混合项，我们可以负担得起 3D 导数的存储。
% 5阶共有 (6*5)/2 = 15 项。
% 3D 存储: 15 * 360*180*40 * 4bytes ~ 3.2 GB。可行。
% 决定：计算全 3D 导数以保证平均态分析的精度。

SA_Mean = gsw_SA_from_SP(Mean_S, P_3D, LON_3D, LAT_3D);
CT_Mean = gsw_CT_from_pt(SA_Mean, Mean_T);

Derivs = cell(MaxOrder, MaxOrder+1);
% 仅计算混合项 (k=1 to n-1)
count = 0;
for n = 2:MaxOrder
    for k = 1:n-1 % k 是 S 的幂次。仅混合项。
        n_T = n - k;
        % 计算 3D 混合导数
        D_val = Engine.calculate_mixed(SA_Mean, CT_Mean, P_3D, n_T, k);
        Derivs{n, k+1} = single(D_val);
        count = count + 1;
    end
end
fprintf('   已计算 %d 个混合导数场\n', count);

%% 5. 遍历计算
fprintf('>> [5/6] 遍历计算混合项...\n');
TaskList = struct('Year', [], 'Month', [], 'TimeVal', [], 'FullPath', []);
task_idx = 0;
for y = Years
    Folder = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = 1:12
        task_idx = task_idx + 1;
        Pattern = sprintf('*%d%02d*.nc', y, m);
        FItem = dir(fullfile(Folder, Pattern));
        TaskList(task_idx).Year = y; TaskList(task_idx).Month = m;
        TaskList(task_idx).TimeVal = y + (m-0.5)/12;
        if ~isempty(FItem), TaskList(task_idx).FullPath = fullfile(FItem(1).folder, FItem(1).name);
        else, TaskList(task_idx).FullPath = ''; end
    end
end
TotalSteps = length(TaskList);

% 初始化结果容器
% 命名规则: Cross_T{nT}S{nS}
% 需要的变量:
% 2阶: T1S1
% 3阶: T2S1, T1S2
% 4阶: T3S1, T2S2, T1S3
% 5阶: T4S1, T3S2, T2S3, T1S4
Variables = {};
for n = 2:MaxOrder
    for k = 1:n-1
        VarName = sprintf('Cross_T%dS%d', n-k, k);
        Variables{end+1} = VarName;
        eval([VarName ' = zeros(Nx, Ny, TotalSteps, ''single'');']);
    end
end

time_vec = [TaskList.TimeVal];
InvFact = zeros(21,1); for i=0:20, InvFact(i+1)=1/factorial(i); end

tic;
for idx = 1:TotalSteps
    FullPath = TaskList(idx).FullPath;
    if isempty(FullPath), continue; end
    
    try
        T = double(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
        S = double(ncread(FullPath, 'salinity',    [1 1 1 1], [Inf Inf Nz 1]));
        MaskNaN = isnan(T) | isnan(S);
        if nanmean(T(:)) > 100, T = T - 273.15; end
        
        SA = gsw_SA_from_SP(S, P_3D, LON_3D, LAT_3D);
        CT = gsw_CT_from_pt(SA, T);
        
        % 扰动量 (相对于平均态)
        dS = SA - SA_Mean; dT = CT - CT_Mean;
        dS(MaskNaN) = 0; dT(MaskNaN) = 0;
        
        % 预计算幂次
        S_Pow = zeros(Nx, Ny, Nz, MaxOrder+1); T_Pow = zeros(Nx, Ny, Nz, MaxOrder+1);
        S_Pow(:,:,:,1) = 1; T_Pow(:,:,:,1) = 1;
        for p = 1:MaxOrder
            S_Pow(:,:,:,p+1) = S_Pow(:,:,:,p) .* dS;
            T_Pow(:,:,:,p+1) = T_Pow(:,:,:,p) .* dT;
        end
        
        % 循环计算混合项
        for n = 2:MaxOrder
            for k = 1:n-1
                nT = n - k; nS = k;
                
                C = InvFact(nS+1) * InvFact(nT+1);
                term = C * Derivs{n, k+1} .* S_Pow(:,:,:,nS+1) .* T_Pow(:,:,:,nT+1);
                
                val = -(sum(term .* dz_perm, 3) / rho0) * 1000;
                
                mask_2d = any(MaskNaN, 3);
                val(mask_2d) = NaN;
                
                VarName = sprintf('Cross_T%dS%d', nT, nS);
                eval([VarName '(:,:,idx) = val;']);
            end
        end
        
        if mod(idx, 10) == 0
            fprintf('   Progress: %.1f%% (%d/%d) | Remain: %.1f min\n', ...
                idx/TotalSteps*100, idx, TotalSteps, (toc/idx)*(TotalSteps-idx)/60);
        end
    catch ME
        fprintf('Error: %s\n', ME.message);
    end
end

%% 6. 保存
fprintf('>> [6/6] 保存结果...\n');
SaveVars = {'Lon', 'Lat', 'time_vec'};
SaveVars = [SaveVars, Variables];
save(SaveName, SaveVars{:}, '-v7.3');
fprintf('>> Done: %s\n', SaveName);
