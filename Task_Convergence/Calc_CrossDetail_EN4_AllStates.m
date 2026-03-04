%% Calc_CrossDetail_EN4_AllStates.m
% =========================================================================
% 功能：计算 EN4 在【平均态】和【标准态】下的 2-8 阶混合详细分布 (28项 CrossDetail)
% 输出：
%   - D:\work\EN4_TSLA_Terms\EN4_CrossDetail_Avg.mat
%   - D:\work\EN4_TSLA_Terms\EN4_CrossDetail_Std.mat
% 引擎：TEOS10_HighOrder_Engine
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
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

%% 2. 初始化引擎 (支持前向求导引擎)
fprintf('>> [1/6] 初始化 TEOS10 高阶引擎 (MaxOrder=%d)...\n', MaxOrder);
CacheFile = fullfile(OutputDir, 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

%% 3. 初始化网格和坐标
fprintf('>> [2/6] 提取网格并建立安全积分空间...\n');
FList_Init = dir(fullfile(DataDir, 'EN.4.2.2.analyses.c13.2005', '*.nc'));
if isempty(FList_Init), error('找不到 EN4 数据文件'); end
SampleFile = fullfile(FList_Init(1).folder, FList_Init(1).name);

Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Level = single(ncread(SampleFile, 'depth'));

Lon = mod(Lon, 360);
DepthIdx = find(Level <= MaxDepth);
CalcDepth = Level(DepthIdx);
Nz = length(CalcDepth);

[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = double(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Nx, Ny, ~] = size(P_3D);

LON_GSW = max(0.5, min(359.5, double(LON_3D)));
LAT_GSW = max(-86, min(90, double(LAT_3D)));
P_GSW   = max(0.1, min(6131, P_3D));

dz = zeros(Nz, 1, 'double');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = 0.5 * (CalcDepth(Nz) - CalcDepth(Nz-1)); 
dz_3D = reshape(dz, 1, 1, Nz); 

%% 4. 生成 28 个变量名称
VarList = {};
for n = 2:MaxOrder
    for k = 1:n-1
        VarList{end+1} = sprintf('Cross_T%dS%d', n-k, k);
    end
end
NumVars = length(VarList);

% 阶乘预处理
InvFact = zeros(21,1,'double'); 
for i=0:20, InvFact(i+1) = 1/factorial(i); end

%% 构建任务列表
fprintf('>> [4/6] 构建文件路径任务调度...\n');
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
time_vec = [TaskList.TimeVal];

%% 5. 定义两种态并计算
States = {'Average', 'StdRef'};
StateSaveNames = {'EN4_CrossDetail_Avg.mat', 'EN4_CrossDetail_Std.mat'};

for st = 1:length(States)
    cur_state = States{st};
    fprintf('\n========== 开始处理 [%s] 态 ==========\n', cur_state);
    
    SavePath = fullfile(OutputDir, StateSaveNames{st});
    if exist(SavePath, 'file')
        fprintf('  --> 已存在，跳过: %s\n', SavePath);
        continue;
    end
    
    % --- Step A: 获得该状态下所有节点的 参考 SA 和 CT ---
    if strcmp(cur_state, 'Average')
        MeanFile = fullfile('D:\work\EN4_TSLA_Terms', 'EN4_Mean_State.mat');
        if ~exist(MeanFile, 'file'), error('缺少平均态文件，请先运行 EN4_4in1.m 生成'); end
        M = load(MeanFile, 'Mean_T', 'Mean_S');
        % 计算 MaskNaN_Mean
        MaskNaN_Mean = isnan(M.Mean_T) | isnan(M.Mean_S);
        
        SA_Ref_3D = double(gsw_SA_from_SP(double(M.Mean_S), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_pt(SA_Ref_3D, double(M.Mean_T))); % EN4 Mean_T 已经是位温了(若原始是的话，按照习惯计算即可)
    else
        % StdRef
        MaskNaN_Mean = false(Nx, Ny, Nz);
        
        SA_Ref_3D = double(gsw_SA_from_SP(Std_Sp * ones(Nx, Ny, Nz, 'double'), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_pt(SA_Ref_3D, Std_Th * ones(Nx, Ny, Nz, 'double')));
    end
    
    % --- Step B: 预计算 28 项的高阶导数空间分布 (Derivs矩阵) ---
    fprintf('>> 计算该态基准背景下的 28 项混合导数矩阵...\n');
    Derivs_Map = containers.Map();
    for n = 2:MaxOrder
        for k = 1:n-1
            nT = n - k; nS = k;
            v_name = sprintf('Cross_T%dS%d', nT, nS);
            Derivs_Map(v_name) = double(Engine.calculate_mixed(SA_Ref_3D, CT_Ref_3D, P_GSW, nT, nS));
        end
    end
    
    % --- Step C: 初始化 28 个时序容器 ---
    DataCubes = struct();
    for v = 1:NumVars
        v_name = VarList{v};
        DataCubes.(v_name) = zeros(Nx, Ny, TotalSteps, 'single');
    end
    
    % --- Step D: 开始循环加载处理每个月的数据 ---
    fprintf('>> 循环计算逐月异常积分分布...\n');
    tic;
    warning('off', 'all');
    for idx = 1:TotalSteps
        FullPath = TaskList(idx).FullPath;
        
        if isempty(FullPath)
            for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
            continue;
        end
        
        try
            T_inst = double(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
            S_inst = double(ncread(FullPath, 'salinity',    [1 1 1 1], [Inf Inf Nz 1]));
            
            % 处理单位与缺测
            if nanmean(T_inst(:)) > 100, T_inst = T_inst - 273.15; end
            
            MaskInvalid = (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0) | isnan(T_inst) | isnan(S_inst);
            MaskValid = ~(MaskNaN_Mean | MaskInvalid);
            idx_valid = find(MaskValid);
            
            if isempty(idx_valid)
                for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
                continue;
            end
            
            DynamicGap = MaskInvalid & ~MaskNaN_Mean;
            HasGap_2D = squeeze(any(DynamicGap, 3)); 
            MaskLand = squeeze(sum(~MaskNaN_Mean, 3)) < 1;
            
            SA_1D = gsw_SA_from_SP(S_inst(idx_valid), P_GSW(idx_valid), LON_GSW(idx_valid), LAT_GSW(idx_valid));
            CT_1D = gsw_CT_from_pt(SA_1D, T_inst(idx_valid));
            
            dTheta_1D = CT_1D - CT_Ref_3D(idx_valid);
            dSA_1D    =  SA_1D - SA_Ref_3D(idx_valid);
            
            for v = 1:NumVars
                v_name = VarList{v};
                parts = sscanf(v_name, 'Cross_T%dS%d');
                nT = parts(1); nS = parts(2);
                
                Combinator = InvFact(nT+1) * InvFact(nS+1);
                Deriv_Valid = Derivs_Map(v_name);
                Deriv_Valid = Deriv_Valid(idx_valid);
                
                C_val_1D = Combinator .* Deriv_Valid .* (dTheta_1D .^ nT) .* (dSA_1D .^ nS);
                
                C_3D = zeros(Nx, Ny, Nz, 'double');
                C_3D(idx_valid) = C_val_1D;
                
                val_C = squeeze(-sum(C_3D .* dz_3D, 3)) / rho0 * 1000; % unit: mm
                val_C(MaskLand | HasGap_2D) = NaN;
                
                DataCubes.(v_name)(:,:,idx) = single(val_C);
            end
        catch ME
            fprintf('   ⚠️ 计算出错 %d: %s\n', idx, ME.message);
            for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
        end
        if mod(idx, 20) == 0, fprintf('   ... 完成 %.1f%%\n', idx/TotalSteps*100); end
    end
    warning('on', 'all');
    
    % --- Step E: 统一坐标命名并保存 ---
    fprintf('>> 保存态结果...\n');
    lon = Lon; lat = Lat; 
    
    MFinal = struct('lon', lon, 'lat', lat, 'time_vec', time_vec);
    for v = 1:NumVars
        v_name = VarList{v};
        MFinal.(v_name) = permute(DataCubes.(v_name), [1, 2, 3]);
    end
    
    save(SavePath, '-struct', 'MFinal', '-v7.3');
    fprintf('🎉 完成并保存至 %s\n', SavePath);
end

fprintf('\n>>> EN4 全部运算完毕！ <<<\n');
