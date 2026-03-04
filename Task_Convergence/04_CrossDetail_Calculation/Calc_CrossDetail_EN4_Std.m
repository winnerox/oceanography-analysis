 %% Calc_CrossDetail_Std.m
% =========================================================================
% 功能：计算【标准态】下的混合偏导数项（逐项分解）
% 阶数：2-8阶 (共 28 项)
% 引擎：TEOS10_HighOrder_Engine
% 参考态：Std_Sp=35, Std_Th=0 (正确定义)
% 
% 输出：EN4_CrossDetail_Std.mat
%   包含变量如: Cross_T1S1 ... Cross_T7S1 等
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
SaveName = fullfile(OutputDir, 'EN4_CrossDetail_Std.mat');
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; 
rho0 = 1035.0; 
Std_Sp = 35.0; 
Std_Th = 0.0;
MaxOrder = 8; 

CacheFile = fullfile(OutputDir, 'TEOS10_Engine_Cache.mat');
fprintf('>> [1/6] 初始化 TEOS10 高阶引擎 (MaxOrder=%d)...\n', MaxOrder);
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false); 

%% 2. 初始化网格
fprintf('>> [2/6] 初始化网格...\n');
FList_Init = dir(fullfile(DataDir, 'EN.4.2.2.analyses.c13.2005', '*.nc'));
if isempty(FList_Init), error('找不到数据文件'); end
SampleFile = fullfile(FList_Init(1).folder, FList_Init(1).name);

Lon = ncread(SampleFile, 'lon'); Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= 2000); CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

%% 3. 预计算标准态场与导数
fprintf('>> [3/6] 预计算标准态场与导数 (Sp=35, Theta=0)...\n');

% 3.1 3D 参考场 (展开中心)
SA_Ref_3D = gsw_SA_from_SP(Std_Sp * ones(Nx, Ny, Nz), P_3D, LON_3D, LAT_3D);
CT_Ref_3D = gsw_CT_from_pt(SA_Ref_3D, Std_Th * ones(Nx, Ny, Nz));

% 3.2 1D 导数系数廓线
P_Prof = gsw_p_from_z(-CalcDepth, 0); 
Nz_prof = length(P_Prof);
SA_Prof = gsw_SA_from_SP(Std_Sp * ones(Nz_prof,1), P_Prof, zeros(Nz_prof,1), zeros(Nz_prof,1)); 
CT_Prof = gsw_CT_from_pt(SA_Prof, Std_Th * ones(Nz_prof,1)); 

Derivs = cell(MaxOrder, MaxOrder+1); 
count = 0;
for n = 2:MaxOrder
    for k = 1:n-1 % k 是 S 的幂次。仅混合项。
        n_T = n - k;
        % 计算原始导数 (快速，不简化)
        D_val = Engine.calculate_mixed(SA_Prof, CT_Prof, P_Prof, n_T, k);
        Derivs{n, k+1} = reshape(D_val, 1, 1, Nz); 
        count = count + 1;
    end
end
fprintf('   已计算 %d 个混合导数项\n', count);

%% 4. 构建任务列表
fprintf('>> [4/6] 任务调度...\n');
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

% 预分配变量
Variables = {};
for n = 2:MaxOrder
    for k = 1:n-1
        VarName = sprintf('Cross_T%dS%d', n-k, k);
        Variables{end+1} = VarName;
        % 动态分配
        eval([VarName ' = zeros(Nx, Ny, TotalSteps, ''single'');']);
    end
end
time_vec = [TaskList.TimeVal];
InvFact = zeros(21,1); for i=0:20, InvFact(i+1)=1/factorial(i); end

%% 5. 遍历计算
fprintf('>> [5/6] 遍历计算 (Total Steps: %d)...\n', TotalSteps);
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
        
        dS = SA - SA_Ref_3D; 
        dT = CT - CT_Ref_3D;
        dS(MaskNaN) = 0; dT(MaskNaN) = 0;
        
        % 预计算幂次
        S_Pow = zeros(Nx, Ny, Nz, MaxOrder+1); T_Pow = zeros(Nx, Ny, Nz, MaxOrder+1);
        S_Pow(:,:,:,1) = 1; T_Pow(:,:,:,1) = 1;
        for p = 1:MaxOrder
            S_Pow(:,:,:,p+1) = S_Pow(:,:,:,p) .* dS;
            T_Pow(:,:,:,p+1) = T_Pow(:,:,:,p) .* dT;
        end
        
        % 循环计算混合项
        mask_2d = any(MaskNaN, 3);
        
        for n = 2:MaxOrder
            for k = 1:n-1
                nS = k; nT = n - k;
                
                C = InvFact(nS+1) * InvFact(nT+1);
                term = C * Derivs{n, k+1} .* S_Pow(:,:,:,nS+1) .* T_Pow(:,:,:,nT+1);
                
                val = -(sum(term .* dz_perm, 3) / rho0) * 1000;
                val(mask_2d) = NaN;
                
                VarName = sprintf('Cross_T%dS%d', nT, nS);
                eval([VarName '(:,:,idx) = val;']);
            end
        end
        
        if mod(idx, 10) == 0
            fprintf('   Progress: %.1f%% | Remain: %.1f min\n', ...
                idx/TotalSteps*100, (toc/idx)*(TotalSteps-idx)/60);
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
