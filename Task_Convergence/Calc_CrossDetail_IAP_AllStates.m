%% Calc_CrossDetail_IAP_AllStates.m
% =========================================================================
% 功能：计算 IAP 在【平均态】和【标准态】下的 2-8 阶混合详细分布 (28项 CrossDetail)
% 输出：
%   - D:\work\IAP_TSLA_Terms\IAP_CrossDetail_Avg.mat
%   - D:\work\IAP_TSLA_Terms\IAP_CrossDetail_Std.mat
% 引擎：原版 TEOS10_General_Engine (保持不动)
% 优化：提取内层循环冗余、预计算指数幂、内存定点复用
% =========================================================================
clear; clc;

%% 1. 配置
fprintf('=== IAP 混合项详细分布计算 ===\n');
fprintf('>> 配置参数...\n');
DataDir_T = 'D:\work\IAP_05_24\TEMP'; 
DataDir_S = 'D:\work\IAP_05_24\SALT'; 
OutputDir = 'D:\work\IAP_mat_data';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024;
Months = 1:12;
TotalSteps = length(Years) * 12;

MaxOrder = 8;
MaxDepth = 2000;
rho0 = double(1035.0);

% 标准态定义
Std_Sp = 35.0;
Std_Th = 0.0;

fprintf('  年份范围: %d-%d\n', min(Years), max(Years));
fprintf('  月份范围: %d-%d\n', min(Months), max(Months));
fprintf('  总时间步数: %d\n', TotalSteps);
fprintf('  最大阶数: %d\n', MaxOrder);
fprintf('  最大深度: %d m\n', MaxDepth);

%% 2. 初始化引擎 (调用你现有的未修改版本)
fprintf('>> [1/6] 初始化 TEOS10 通用引擎 (MaxOrder=%d)...\n', MaxOrder);
Engine = TEOS10_General_Engine(MaxOrder);

%% 3. 初始化网格和坐标
fprintf('>> [2/6] 提取网格并建立安全积分空间...\n');
FList_Init = dir(fullfile(DataDir_T, '*2005*01*.nc')); 
if isempty(FList_Init), error('找不到 IAP 温度数据文件'); end
SampleFile = fullfile(FList_Init(1).folder, FList_Init(1).name);

fprintf('  读取网格信息...\n');
Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Level = single(ncread(SampleFile, 'depth_std')); 

Lon = mod(Lon, 360);
DepthIdx = find(Level <= MaxDepth);
CalcDepth = Level(DepthIdx);
Nz = length(CalcDepth);

[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon);
P_3D = double(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Ny, ~, Nx] = size(P_3D);

LON_GSW = max(0.5, min(359.5, double(LON_3D)));
LAT_GSW = max(-86, min(90, double(LAT_3D)));
P_GSW   = max(0.1, min(6131, P_3D));

dz = zeros(Nz, 1, 'double');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = 0.5 * (CalcDepth(Nz) - CalcDepth(Nz-1)); 
dz_3D = reshape(dz, 1, Nz, 1);

fprintf('  网格信息: %d (lat) x %d (lon) x %d (depth)\n', Ny, Nx, Nz); 

%% 4. 生成 28 个变量名称
fprintf('>> [3/6] 生成混合项变量名称...\n');
VarList = {};
for n = 2:MaxOrder
    for k = 1:n-1
        VarList{end+1} = sprintf('Cross_T%dS%d', n-k, k);
    end
end
NumVars = length(VarList);

fprintf('  生成 %d 个混合项变量\n', NumVars);
for v = 1:NumVars
    if mod(v, 7) == 0
        fprintf('  %s\n', VarList{v});
    else
        fprintf('  %s ', VarList{v});
    end
end

% 阶乘预处理
InvFact = zeros(21,1,'double'); 
for i=0:20, InvFact(i+1) = 1/factorial(i); end
fprintf('  完成阶乘预处理\n');

%% 构建任务列表
fprintf('>> [4/6] 构建文件路径任务调度...\n');
TaskList = struct('Year', [], 'Month', [], 'TimeVal', [], 'FullPath_T', [], 'FullPath_S', []);
task_idx = 0;
valid_tasks = 0;

for y = Years
    for m = Months
        task_idx = task_idx + 1;
        
        Pattern_T = sprintf('*year_%d_month_%02d*.nc', y, m);
        Pattern_S = sprintf('*year_%d_month_%02d*.nc', y, m);
        
        FItem_T = dir(fullfile(DataDir_T, Pattern_T));
        FItem_S = dir(fullfile(DataDir_S, Pattern_S));
        
        TaskList(task_idx).Year = y; TaskList(task_idx).Month = m;
        TaskList(task_idx).TimeVal = y + (m-0.5)/12;
        
        if ~isempty(FItem_T), TaskList(task_idx).FullPath_T = fullfile(FItem_T(1).folder, FItem_T(1).name); else, TaskList(task_idx).FullPath_T = ''; end
        if ~isempty(FItem_S), TaskList(task_idx).FullPath_S = fullfile(FItem_S(1).folder, FItem_S(1).name); else, TaskList(task_idx).FullPath_S = ''; end
        
        if ~isempty(FItem_T) && ~isempty(FItem_S)
            valid_tasks = valid_tasks + 1;
        end
    end
end
time_vec = [TaskList.TimeVal];
fprintf('  完成任务列表构建，有效任务数: %d/%d\n', valid_tasks, TotalSteps);

%% 5. 定义两种态并计算
States = {'Average', 'StdRef'};
StateSaveNames = {'IAP_CrossDetail_Avg.mat', 'IAP_CrossDetail_Std.mat'};

for st = 1:length(States)
    cur_state = States{st};
    fprintf('\n========== 开始处理 [%s] 态 ==========\n', cur_state);
    
    SavePath = fullfile(OutputDir, StateSaveNames{st});
    if exist(SavePath, 'file')
        fprintf('  --> 已存在，跳过: %s\n', SavePath);
        continue;
    end
    
    % --- Step A: 获得该状态下所有节点的 参考 SA 和 CT ---
    fprintf('>> [5/6] 计算参考状态...\n');
    if strcmp(cur_state, 'Average')
        MeanFile = fullfile('D:\work\IAP_mat_data', 'IAP_Mean_State_05_24.mat');
        if ~exist(MeanFile, 'file'), error('缺少平均态文件，请先运行相应的均态生成'); end
        fprintf('  加载平均态文件: %s\n', MeanFile);
        M = load(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean');
        MaskNaN_Mean = M.MaskNaN_Mean;
        
        fprintf('  计算参考 SA 和 CT...\n');
        SA_Ref_3D = double(gsw_SA_from_SP(double(M.Mean_S), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_pt(SA_Ref_3D, double(M.Mean_T))); 
    else
        MaskNaN_Mean = false(Ny, Nz, Nx); 
        fprintf('  计算标准态 SA 和 CT...\n');
        SA_Ref_3D = double(gsw_SA_from_SP(Std_Sp * ones(Ny, Nz, Nx, 'double'), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_pt(SA_Ref_3D, Std_Th * ones(Ny, Nz, Nx, 'double')));
    end
    
    % --- Step B: 调用现成引擎，一次性获取所有导数空间分布 ---
    fprintf('>> 计算该态基准背景下的所有混合导数矩阵...\n');
    tic;
    [~, ~, Derivs_Cross_All] = Engine.calculate_all_derivatives(SA_Ref_3D, CT_Ref_3D, P_GSW, MaxOrder);
    toc;
    
    % --- Step C: 初始化时序容器与预处理参数 (提速核心：前置运算) ---
    fprintf('>> 初始化数据容器和预处理参数...\n');
    DataCubes = struct();
    nT_list = zeros(NumVars, 1);
    nS_list = zeros(NumVars, 1);
    Combinator_list = zeros(NumVars, 1);
    Derivs_Cell = cell(NumVars, 1); 
    
    for v = 1:NumVars
        v_name = VarList{v};
        DataCubes.(v_name) = zeros(Ny, Nx, TotalSteps, 'single');
        
        % 预处理提取 nT 和 nS (干掉循环内的 sscanf)
        parts = sscanf(v_name, 'Cross_T%dS%d');
        nT_list(v) = parts(1); 
        nS_list(v) = parts(2);
        
        % 预计算组合常数
        Combinator_list(v) = InvFact(nT_list(v)+1) * InvFact(nS_list(v)+1);
        
        % 直接将引擎计算好的 3D 矩阵存入 Cell，后续按索引极速读取
        Derivs_Cell{v} = Derivs_Cross_All{nT_list(v), nS_list(v)};
    end
    fprintf('  完成数据容器初始化\n');
    
    % --- Step D: 开始循环加载处理每个月的数据 ---
    fprintf('>> [6/6] 循环计算逐月异常积分分布...\n');
    tic;
    warning('off', 'all');
    
    % [提速核心：内存复用] 预分配 3D 缓存矩阵，避免 28 x TotalSteps 次的内存开辟
    C_3D_buffer = zeros(Ny, Nz, Nx, 'double');
    
    for idx = 1:TotalSteps
        FullPath_T = TaskList(idx).FullPath_T;
        FullPath_S = TaskList(idx).FullPath_S;
        
        if isempty(FullPath_T) || isempty(FullPath_S)
            fprintf('  ⚠️ 跳过时间步 %d: 数据文件缺失\n', idx);
            for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
            continue;
        end
        
        try
            fprintf('  处理时间步 %d/%d: %d年%d月\n', idx, TotalSteps, TaskList(idx).Year, TaskList(idx).Month);
            
            T_inst = double(ncread(FullPath_T, 'temp', [1 1 1], [Nz Inf Inf])); 
            S_inst = double(ncread(FullPath_S, 'salinity', [1 1 1], [Nz Inf Inf]));
            T_inst = permute(T_inst, [3, 1, 2]);
            S_inst = permute(S_inst, [3, 1, 2]);
            
            if nanmean(T_inst(:)) > 100, T_inst = T_inst - 273.15; end
            
            MaskInvalid = (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0) | isnan(T_inst) | isnan(S_inst);
            MaskValid = ~(MaskNaN_Mean | MaskInvalid);
            idx_valid = find(MaskValid);
            
            if isempty(idx_valid)
                fprintf('  ⚠️ 时间步 %d: 无有效数据\n', idx);
                for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
                continue;
            end
            
            DynamicGap = MaskInvalid & ~MaskNaN_Mean;
            HasGap_2D = squeeze(any(DynamicGap, 2)); 
            MaskLand = squeeze(sum(~MaskNaN_Mean, 2)) < 1;
            
            SA_1D = gsw_SA_from_SP(S_inst(idx_valid), P_GSW(idx_valid), LON_GSW(idx_valid), LAT_GSW(idx_valid));
            CT_1D = gsw_CT_from_pt(SA_1D, T_inst(idx_valid));
            
            dTheta_1D = CT_1D - CT_Ref_3D(idx_valid);
            dSA_1D    =  SA_1D - SA_Ref_3D(idx_valid);
            
            % [提速核心：算法降维] 预计算幂次方，取代内部耗时的指数运算 (.^)
            pow_dTheta = ones(length(idx_valid), MaxOrder, 'double');
            pow_dSA    = ones(length(idx_valid), MaxOrder, 'double');
            pow_dTheta(:, 1) = dTheta_1D;
            pow_dSA(:, 1)    = dSA_1D;
            for p = 2:MaxOrder
                pow_dTheta(:, p) = pow_dTheta(:, p-1) .* dTheta_1D;
                pow_dSA(:, p)    = pow_dSA(:, p-1) .* dSA_1D;
            end
            
            % 最内层循环现在只有纯净的标量乘法和数组赋值
            for v = 1:NumVars
                v_name = VarList{v};
                nT = nT_list(v); 
                nS = nS_list(v);
                
                % 1. 提取预计算好的底层导数参数
                Deriv_Valid = Derivs_Cell{v}(idx_valid);
                
                % 2. 纯乘法运算代替复杂的幂运算和除法
                C_val_1D = Combinator_list(v) .* Deriv_Valid .* pow_dTheta(:, nT) .* pow_dSA(:, nS);
                
                % 3. 复用预分配内存映射回 3D 空间
                C_3D_buffer(idx_valid) = C_val_1D;
                
                % 4. 垂直方向快速积分
                val_C = squeeze(-sum(C_3D_buffer .* dz_3D, 2)) * (1000 / rho0); 
                val_C(MaskLand | HasGap_2D) = NaN;
                
                DataCubes.(v_name)(:,:,idx) = single(val_C);
                
                % 5. 快速定点清空缓冲区
                C_3D_buffer(idx_valid) = 0;
            end
        catch ME
            fprintf('   ⚠️ 计算出错 %d: %s\n', idx, ME.message);
            for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
        end
        if mod(idx, 5) == 0, fprintf('   ... 完成 %.1f%%\n', idx/TotalSteps*100); end
    end
    warning('on', 'all');
    toc;
    
    % --- Step E: 统一坐标命名并保存 ---
    fprintf('>> 保存态结果...\n');
    lon = Lon; lat = Lat; 
    
    MFinal = struct('lon', lon, 'lat', lat, 'time_vec', time_vec);
    for v = 1:NumVars
        v_name = VarList{v};
        MFinal.(v_name) = permute(DataCubes.(v_name), [2, 1, 3]);
    end
    
    save(SavePath, '-struct', 'MFinal', '-v7.3');
    fprintf('🎉 完成并保存至 %s\n', SavePath);
end

fprintf('\n>>> IAP 全部运算完毕！ <<<\n');