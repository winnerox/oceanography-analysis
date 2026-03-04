%% Calc_Cross_Terms_EN4_Std.m
% =========================================================================
% 功能：计算 EN4【标准态】下的混合项 (Cross Terms, 2-8阶)
% 输出：D:\work\MAT_Data\EN4_Cross_Terms_1to8_StdRef.mat
% 极致优化：
%   1. 一维展平跳过陆地：仅对有效水体计算TEOS-10
%   2. 避免4D内存杀手：不预计算S_Pow/T_Pow，现场计算幂次
%   3. 强制单精度：所有变量显式single类型
%   4. 一阶冗余：从n=2开始
% =========================================================================
clear; clc;
addpath('D:\work');

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\MAT_Data';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'EN4_Cross_Terms_1to8_StdRef.mat');

Years = 2005:2024; 
rho0 = single(1035.0); 
Std_Sp = single(35.0);
Std_Th = single(0.0);
MaxOrder = 8;
MaxDepth = 2000;

%% 2. 初始化 TEOS10 高阶引擎
CacheFile = fullfile(OutputDir, 'TEOS10_Engine_Cache.mat');
fprintf('[1/6] 初始化 TEOS10 高阶引擎 (MaxOrder=%d)...\n', MaxOrder);
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

%% 3. 初始化网格
fprintf('[2/6] 初始化网格...\n');
FList_Init = dir(fullfile(DataDir, 'EN.4.2.2.analyses.c13.2005', '*.nc'));
if isempty(FList_Init), error('找不到数据文件'); end
SampleFile = fullfile(FList_Init(1).folder, FList_Init(1).name);

Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, 'depth')); 
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

%% 4. 预计算标准态场与导数（使用3D场，保持空间一致性）
fprintf('[3/6] 预计算标准态场与导数 (Sp=35, Theta=0)...\n');

SA_Ref_3D = single(gsw_SA_from_SP(double(Std_Sp * ones(Nx, Ny, Nz, 'single')), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref_3D = single(gsw_CT_from_pt(SA_Ref_3D, double(Std_Th * ones(Nx, Ny, Nz, 'single'))));

Derivs = cell(MaxOrder, MaxOrder+1);
for n = 1:MaxOrder
    for k = 0:n
        n_S = k; n_T = n - k;
        D_val = Engine.calculate_mixed(SA_Ref_3D, CT_Ref_3D, P_3D, n_T, n_S);
        Derivs{n, k+1} = single(D_val);
    end
end

InvFact = zeros(21,1,'single'); 
for i=0:20, InvFact(i+1) = single(1/factorial(i)); end

%% 5. 构建任务列表
fprintf('[4/6] 构建文件索引列表...\n');
TaskList = struct('Year', [], 'Month', [], 'TimeVal', [], 'FullPath', []);
task_idx = 0;

for y = Years
    Folder = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = 1:12
        task_idx = task_idx + 1;
        Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        FItem = dir(fullfile(Folder, Pattern));
        
        TaskList(task_idx).Year = y;
        TaskList(task_idx).Month = m;
        TaskList(task_idx).TimeVal = y + (m-0.5)/12;
        
        if ~isempty(FItem)
            TaskList(task_idx).FullPath = fullfile(FItem(1).folder, FItem(1).name);
        else
            TaskList(task_idx).FullPath = '';
        end
    end
end
TotalSteps = length(TaskList);

%% 6. 遍历计算
fprintf('[5/6] 开始遍历计算 (Total Steps: %d)...\n', TotalSteps);

Cross_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_vec = [TaskList.TimeVal];

tic;
for idx = 1:TotalSteps
    FullPath = TaskList(idx).FullPath;
    if isempty(FullPath), continue; end
    
    try
        T = single(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
        S = single(ncread(FullPath, 'salinity',    [1 1 1 1], [Inf Inf Nz 1]));
        
        if mean(T(:), 'omitnan') > 100, T = T - single(273.15); end
        MaskNaN = isnan(T) | isnan(S) | (T <= -100) | (T >= 100) | (S < 0) | (S >= 100);
        
        SA = single(gsw_SA_from_SP(double(S), double(P_3D), double(LON_3D), double(LAT_3D)));
        CT = single(gsw_CT_from_pt(SA, double(T)));
        
        dS = SA - SA_Ref_3D;
        dT = CT - CT_Ref_3D;
        dS(MaskNaN) = 0;
        dT(MaskNaN) = 0;
        
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
        
        if mod(idx, 10) == 0
            avg_time = toc / idx;
            remain_time = avg_time * (TotalSteps - idx);
            fprintf('   进度: %.1f%% (%d/%d) | 剩余: %.1f min\n', ...
                idx/TotalSteps*100, idx, TotalSteps, remain_time/60);
        end
        
    catch ME
        fprintf('Error in %d-%d: %s\n', TaskList(idx).Year, TaskList(idx).Month, ME.message);
    end
end

%% 7. 保存结果
fprintf('[6/6] 保存结果...\n');
lon = Lon; lat = Lat;
save(SaveName, 'Cross_AllOrders', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('完成! 文件: %s\n', SaveName);
