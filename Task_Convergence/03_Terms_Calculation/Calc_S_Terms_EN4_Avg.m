%% Calc_S_Terms_EN4_Avg.m
% =========================================================================
% 功能：计算 EN4【平均态】下的纯盐度项 (Pure S Terms, 1-8阶)
% 输出：D:\work\MAT_Data\EN4_S_Terms_1to8_Average.mat
% 极致优化：
%   1. 掩码逻辑：只对有效深度求和，保留浅海数据
%   2. 强制单精度：所有变量显式single类型
% =========================================================================
clear; clc;
addpath('D:\work');

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\MAT_Data';
MeanStateFile = 'D:\work\EN4_TSLA_Terms\EN4_Mean_State.mat';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
SaveName = fullfile(OutputDir, 'EN4_S_Terms_1to8_Average.mat');

Years = 2005:2024; 
rho0 = single(1035.0);       
MaxDepth = 2000;     
MaxOrder = 8;

%% 2. 加载平均态
fprintf('[1/6] 加载平均态数据...\n');
if ~exist(MeanStateFile, 'file'), error('找不到平均态文件'); end
load(MeanStateFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');
Mean_T = single(Mean_T); Mean_S = single(Mean_S);
Lon = single(Lon); Lat = single(Lat);
CalcDepth = single(Depth(DepthIdx));

%% 3. 初始化网格
fprintf('[2/6] 初始化网格...\n');
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

%% 4. 初始化 TEOS10 引擎
fprintf('[3/6] 初始化 TEOS-10 引擎...\n');
CacheFile = fullfile(OutputDir, 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

SA_Ref = single(gsw_SA_from_SP(double(Mean_S), double(P_3D), double(LON_3D), double(LAT_3D)));
CT_Ref = single(gsw_CT_from_pt(SA_Ref, double(Mean_T)));

Derivs_S = cell(MaxOrder, 1);
for n = 1:MaxOrder
    D_val = Engine.calculate_mixed(SA_Ref, CT_Ref, P_3D, 0, n);
    Derivs_S{n} = single(D_val);
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

SSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
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
        
        dS = SA - SA_Ref;
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
save(SaveName, 'SSLA_AllOrders', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('完成! 文件: %s\n', SaveName);
