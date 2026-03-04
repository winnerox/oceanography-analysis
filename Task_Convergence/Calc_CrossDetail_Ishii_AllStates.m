%% Calc_CrossDetail_Ishii_AllStates.m
% =========================================================================
% 功能：计算 Ishii 在【平均态】和【标准态】下的 2-8 阶混合详细分布 (28项 CrossDetail)
% 输出：
%   - D:\work\Ishii_mat_data\Ishii_CrossDetail_Avg.mat
%   - D:\work\Ishii_mat_data\Ishii_CrossDetail_Std.mat
% 引擎：TEOS10_HighOrder_Engine
% 参考态：
%   - 平均态: 使用 2005-2024 本地网格多年平均 T/S 场
%   - 标准态: 使用 全局统一常数 Sp=35, Theta=0
%
% 注意：本脚本采用了 Ishii_4in1 中经受过极速优化与防相消截断误差处理的算法结构！
% =========================================================================
clear; clc;

%% 1. 配置
DataRoot = 'D:\work\Ishii_05_24';
TempDir = fullfile(DataRoot, 'Temperature');
SaltDir = fullfile(DataRoot, 'Salinity');
OutputDir = 'D:\work\Ishii_mat_data'; 
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
sample_file = fullfile(TempDir, 'temp.2005.nc');
Lon = single(ncread(sample_file, 'longitude'));
Lat = single(ncread(sample_file, 'latitude'));
Level = single(ncread(sample_file, 'level'));

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

%% 4. 生成 28 个变量名称
VarList = {};
for n = 2:MaxOrder
    for k = 1:n-1
        VarList{end+1} = sprintf('Cross_T%dS%d', n-k, k);
    end
end
NumVars = length(VarList);

% 预构建时间轴
Time_Axis = zeros(TotalSteps, 1);
for y_idx = 1:length(Years)
    for m = 1:12
        idx = (y_idx - 1) * 12 + m;
        Time_Axis(idx) = Years(y_idx) + (m - 0.5) / 12;
    end
end

% 阶乘预处理
InvFact = zeros(21,1,'double'); 
for i=0:20, InvFact(i+1) = 1/factorial(i); end

%% 5. 定义两种态并计算
States = {'Average', 'StdRef'};
StateSaveNames = {'Ishii_CrossDetail_Avg.mat', 'Ishii_CrossDetail_Std.mat'};

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
        MeanFile = fullfile(OutputDir, 'Ishii_Mean_State_05_24.mat');
        if ~exist(MeanFile, 'file'), error('缺少平均态文件，请先运行 Ishii_4in1.m'); end
        M = load(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean');
        MaskNaN_Mean = M.MaskNaN_Mean;
        
        SA_Ref_3D = double(gsw_SA_from_SP(double(M.Mean_S), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_t(SA_Ref_3D, double(M.Mean_T), P_GSW));
    else
        % StdRef
        % 对于标准态其实不需要 MaskNaN_Mean 在预计算时应用，因为在处理每时每刻的时会动态丢弃
        MaskNaN_Mean = false(Ny, Nz, Nx);
        
        SA_Ref_3D = double(gsw_SA_from_SP(Std_Sp * ones(Ny, Nz, Nx, 'double'), P_GSW, LON_GSW, LAT_GSW));
        CT_Ref_3D = double(gsw_CT_from_pt(SA_Ref_3D, Std_Th * ones(Ny, Nz, Nx, 'double')));
    end
    
    % --- Step B: 预计算 28 项的高阶导数空间分布 (Derivs矩阵) ---
    fprintf('>> 计算该态基准背景下的 28 项混合导数矩阵...\n');
    Derivs_Map = containers.Map();
    for n = 2:MaxOrder
        for k = 1:n-1
            nT = n - k; nS = k;
            v_name = sprintf('Cross_T%dS%d', nT, nS);
            % 这里使用 3D 矩阵直接算会很大，但利用 GSW 高阶底层可以直接算出所有格点在这个参考态的偏导率
            Derivs_Map(v_name) = double(Engine.calculate_mixed(SA_Ref_3D, CT_Ref_3D, P_GSW, nT, nS));
        end
    end
    
    % --- Step C: 初始化 28 个时序容器 ---
    DataCubes = struct();
    for v = 1:NumVars
        v_name = VarList{v};
        DataCubes.(v_name) = zeros(Ny, Nx, TotalSteps, 'single');
    end
    
    % --- Step D: 开始循环加载处理每个月的数据 ---
    fprintf('>> 循环计算逐月异常积分分布...\n');
    tic;
    warning('off', 'all');
    for y = Years
        FileT = fullfile(TempDir, sprintf('temp.%d.nc', y));
        FileS = fullfile(SaltDir, sprintf('sal.%d.nc', y));
        
        if ~exist(FileT, 'file') || ~exist(FileS, 'file')
            for m = 1:12
                idx = (y - Years(1)) * 12 + m; 
                for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
            end
            continue;
        end
        
        try
            T_inst_full = double(permute(ncread(FileT, 'temp', [1 1 1 1], [Inf Inf Nz Inf]), [2, 3, 1, 4]));
            S_inst_full = double(permute(ncread(FileS, 'sal', [1 1 1 1], [Inf Inf Nz Inf]), [2, 3, 1, 4]));
            
            for m = 1:12
                idx = (y - Years(1)) * 12 + m; 
                
                T_inst = T_inst_full(:, :, :, m);
                S_inst = S_inst_full(:, :, :, m);
                
                MaskInvalid = (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0) | isnan(T_inst) | isnan(S_inst);
                MaskValid = ~(MaskNaN_Mean | MaskInvalid);
                idx_valid = find(MaskValid);
                
                if isempty(idx_valid)
                    for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
                    continue;
                end
                
                DynamicGap = MaskInvalid & ~MaskNaN_Mean;
                HasGap_2D = squeeze(any(DynamicGap, 2)); 
                MaskLand = squeeze(sum(~MaskNaN_Mean, 2)) < 1;
                
                % 本月状态
                SA_1D = gsw_SA_from_SP(S_inst(idx_valid), P_GSW(idx_valid), LON_GSW(idx_valid), LAT_GSW(idx_valid));
                CT_1D = gsw_CT_from_t(SA_1D, T_inst(idx_valid), P_GSW(idx_valid));
                
                % 与参考态的差
                dTheta_1D = CT_1D - CT_Ref_3D(idx_valid);
                dSA_1D    =  SA_1D - SA_Ref_3D(idx_valid);
                
                % 对 28 个变量分别积分
                for v = 1:NumVars
                    v_name = VarList{v};
                    parts = sscanf(v_name, 'Cross_T%dS%d');
                    nT = parts(1); nS = parts(2);
                    
                    Combinator = InvFact(nT+1) * InvFact(nS+1);
                    Deriv_Valid = Derivs_Map(v_name);
                    Deriv_Valid = Deriv_Valid(idx_valid);
                    
                    C_val_1D = Combinator .* Deriv_Valid .* (dTheta_1D .^ nT) .* (dSA_1D .^ nS);
                    
                    C_3D = zeros(Ny, Nz, Nx, 'double');
                    C_3D(idx_valid) = C_val_1D;
                    
                    % 深度积分
                    val_C = squeeze(-sum(C_3D .* dz_3D, 2)) / rho0 * 1000; % unit: mm
                    val_C(MaskLand | HasGap_2D) = NaN;
                    
                    DataCubes.(v_name)(:,:,idx) = single(val_C);
                end
            end
        catch ME
            fprintf('   ⚠️ 计算出错 %d 年: %s\n', y, ME.message);
            for m = 1:12
                idx = (y - Years(1)) * 12 + m; 
                for v = 1:NumVars, DataCubes.(VarList{v})(:,:,idx) = NaN; end
            end
        end
        fprintf('   ... 完成 %d 年解析\n', y);
    end
    warning('on', 'all');
    
    % --- Step E: 统一坐标命名并保存 ---
    fprintf('>> 保存态结果...\n');
    time_vec = Time_Axis;
    lon = Lon; lat = Lat; 
    
    MFinal = struct('lon', lon, 'lat', lat, 'time_vec', time_vec);
    for v = 1:NumVars
        v_name = VarList{v};
        % 根据画图习惯重排为 [lon, lat, time] 或原先的 [lat, lon, time]
        % 这里原数据空间是 [lat, lon]，为了和已有保持一致我们转一下
        MFinal.(v_name) = permute(DataCubes.(v_name), [2, 1, 3]);
    end
    
    save(SavePath, '-struct', 'MFinal', '-v7.3');
    fprintf('🎉 完成并保存至 %s\n', SavePath);
end

fprintf('\n>>> 全部运算完毕！ <<<\n');
