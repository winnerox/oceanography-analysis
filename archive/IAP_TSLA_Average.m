%% IAP_TSLA_Average.m
% =========================================================================
% IAP (CAS) 海洋数据 TSLA 计算脚本 [平均态版本]
% =========================================================================
% 功能说明：
%   基于 IAP (中国科学院大气物理研究所) 的温度/盐度格点数据，
%   计算热比容海平面异常 (Thermosteric Sea Level Anomaly, TSLA) 的各阶贡献。
%
% 核心算法：
%   1. Phase 1: 计算 2005-2024 年的多波束长期平均态 (Mean State)。
%   2. Phase 2: 基于平均态，利用 TEOS-10 计算海水状态方程的高阶导数。
%   3. Phase 3: 计算每月的 TSLA 及其泰勒展开项 (1-8 阶)。
%
% 输入数据：
%   - 温度: IAPv4_Temp_monthly_... (单位: °C)
%   - 盐度: IAPv2_Salinity_monthly_... (单位: PSU)
%
% 输出数据：
%   - IAP_TSLA_Terms_1to8_Average.mat
% =========================================================================
clear; clc;

%% ========================================================================
%% 1. 基础参数与路径配置
%% ========================================================================
DataRoot = 'D:\work\IAP_05_24';
TempDir = fullfile(DataRoot, 'TEMP');
SaltDir = fullfile(DataRoot, 'SALT');
OutputDir = 'D:\work\MAT_Data';

% 时间范围
Years = 2005:2024;
Months = 1:12;

% 计算参数
MaxOrder = 8;        % 泰勒展开最高阶数
MaxDepth = 2000;     % 垂直积分深度 (米)
rho0 = 1035.0;       % 参考密度 (kg/m^3)

if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% ========================================================================
%% 2. 初始化核心引擎
%% ========================================================================
fprintf('>> [1/6] 初始化 TEOS-10 热力学引擎...\n');
try
    Engine = TEOS10_General_Engine();
catch
    error('❌ 无法加载 TEOS10_General_Engine，请检查相关脚本是否在路径中。');
end

%% ========================================================================
%% 3. 自动侦测网格信息
%% ========================================================================
fprintf('>> [2/6] 正在侦测网格与变量定义...\n');

% 3.1 搜索首个有效文件作为样本
SampleFile = '';
for y = Years
    for m = 1:12
        Pattern = sprintf('*Temp*year_%d_month_%02d.nc', y, m);
        FileList = dir(fullfile(TempDir, Pattern));
        if ~isempty(FileList)
            SampleFile = fullfile(FileList(1).folder, FileList(1).name);
            break;
        end
    end
    if ~isempty(SampleFile), break; end
end

if isempty(SampleFile)
    error('❌ 未找到任何 IAP 温度文件，请检查路径: %s', TempDir);
end
fprintf('   样本文件: %s\n', SampleFile);

% 3.2 读取元数据并映射变量名
% 屏蔽 NetCDF 版本警告
warning('off', 'all');
Info = ncinfo(SampleFile);
warning('on', 'all');

VarNames = {Info.Variables.Name};
fprintf('   检测变量: %s\n', strjoin(VarNames, ', '));

% 智能匹配变量名
if any(strcmpi('temp', VarNames)), TVar = 'temp';
elseif any(strcmpi('temperature', VarNames)), TVar = 'temperature';
else, error('❌ 无法识别温度变量名'); end

if any(strcmpi('lon', VarNames)), LonVar = 'lon';
elseif any(strcmpi('longitude', VarNames)), LonVar = 'longitude';
else, error('❌ 无法识别经度变量名'); end

if any(strcmpi('lat', VarNames)), LatVar = 'lat';
elseif any(strcmpi('latitude', VarNames)), LatVar = 'latitude';
else, error('❌ 无法识别纬度变量名'); end

if any(strcmpi('depth_std', VarNames)), DepthVar = 'depth_std';
elseif any(strcmpi('depth', VarNames)), DepthVar = 'depth';
elseif any(strcmpi('lev', VarNames)), DepthVar = 'lev';
else, error('❌ 无法识别深度变量名'); end

fprintf('   变量映射: T=[%s], Lon=[%s], Lat=[%s], Depth=[%s]\n', TVar, LonVar, LatVar, DepthVar);

% 3.3 读取基础坐标
Lon = double(ncread(SampleFile, LonVar));
Lat = double(ncread(SampleFile, LatVar));
Depth = double(ncread(SampleFile, DepthVar));

% 确定盐度变量名 (假设文件结构对称)
SampleS_File = strrep(SampleFile, 'TEMP', 'SALT');
SampleS_File = strrep(SampleS_File, 'Temp', 'Salinity');
SampleS_File = strrep(SampleS_File, 'IAPv4', 'IAPv2'); % IAP 版本号差异修正
InfoS = ncinfo(SampleS_File);
VarNamesS = {InfoS.Variables.Name};

if any(strcmpi('salinity', VarNamesS)), SVar = 'salinity';
elseif any(strcmpi('sal', VarNamesS)), SVar = 'sal';
elseif any(strcmpi('salt', VarNamesS)), SVar = 'salt';
else, error('❌ 无法识别盐度变量名'); end
fprintf('   盐度变量: [%s]\n', SVar);

% 3.4 构建 3D 计算网格
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
Nz = length(CalcDepth);

fprintf('   计算网格: %d(Lon) x %d(Lat) x %d(Depth)\n', length(Lon), length(Lat), Nz);

% 使用 ndgrid 生成 3D 坐标矩阵 (注意 MATLAB 维度顺序)
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);

% 【关键】GSW SAAR 参考数据纬度范围是 [-86, 90]，超出范围会导致索引错误
LON_3D = max(0.5, min(359.5, LON_3D));
LAT_3D = max(-86, min(90, LAT_3D));

% 计算静水压力 (用于 GSW 函数)
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
P_3D = max(0.1, min(6131, P_3D));  % GSW 压力范围

% 预计算垂直积分权重 (梯形法则)
dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_3D = reshape(dz, 1, 1, []); % 维度适配 [1, 1, Nz]

%% ========================================================================
%% 4. Phase 1: 计算多年平均态 (Mean State)
%% ========================================================================
fprintf('>> [3/6] Phase 1 - 计算 20 年平均气候态...\n');

Sum_T = zeros(size(P_3D));
Sum_S = zeros(size(P_3D));
Count = zeros(size(P_3D), 'uint16'); % 计数器用于无效值剔除

warning('off', 'all');
for y = Years
    for m = Months
        % 构造文件名
        PatT = sprintf('*Temp*year_%d_month_%02d.nc', y, m);
        PatS = sprintf('*Salinity*year_%d_month_%02d.nc', y, m);
        
        dT = dir(fullfile(TempDir, PatT));
        dS = dir(fullfile(SaltDir, PatS));
        
        if isempty(dT) || isempty(dS), continue; end
        
        FileT = fullfile(dT(1).folder, dT(1).name);
        FileS = fullfile(dS(1).folder, dS(1).name);
        
        try
            T_raw = double(ncread(FileT, TVar));
            S_raw = double(ncread(FileS, SVar));
            
            % [维度修正] IAP 原始数据可能是 [Depth, Lon, Lat] 或其他顺序
            % 这里的逻辑是确保转换为 [Lon, Lat, Depth]
            if size(T_raw, 1) < size(T_raw, 2) && size(T_raw, 1) < size(T_raw, 3)
                T_raw = permute(T_raw, [2, 3, 1]);
                S_raw = permute(S_raw, [2, 3, 1]);
            end
            
            % 深度截取
            if size(T_raw, 3) > Nz
                T_raw = T_raw(:,:,1:Nz);
                S_raw = S_raw(:,:,1:Nz);
            end
            
            % 异常值过滤 (IAP 将 999 设为无效值)
            T_raw(T_raw > 100 | T_raw < -10) = NaN;
            S_raw(S_raw > 100 | S_raw < 0) = NaN;
            
            % 累加
            valid = ~isnan(T_raw) & ~isnan(S_raw);
            Sum_T(valid) = Sum_T(valid) + T_raw(valid);
            Sum_S(valid) = Sum_S(valid) + S_raw(valid);
            Count(valid) = Count(valid) + 1;
        catch
            continue;
        end
    end
    fprintf('   ... 年份 %d 平均态累加完成\n', y);
end
warning('on', 'all');

% 计算平均值
Mean_T = Sum_T ./ double(Count);
Mean_S = Sum_S ./ double(Count);
Mean_T(Count == 0) = NaN; % 无效点置空
Mean_S(Count == 0) = NaN;

clear Sum_T Sum_S Count; % 清理内存

%% ========================================================================
%% 5. Phase 2: 计算背景场热力学导数
%% ========================================================================
fprintf('>> [4/6] Phase 2 - 计算背景场高阶导数 (1-%d 阶)...\n', MaxOrder);

valid_grid = ~isnan(Mean_T) & ~isnan(Mean_S);
fprintf('   有效计算格点数: %d\n', sum(valid_grid(:)));

% 提取有效点进行向量化计算
S_vec = double(Mean_S(valid_grid));
T_vec = double(Mean_T(valid_grid));
P_vec = double(P_3D(valid_grid));
Lon_vec = double(LON_3D(valid_grid));
Lat_vec = double(LAT_3D(valid_grid));

% [坐标修正] GSW SAAR 参考数据纬度范围是 [-86, 90]
Lon_vec = mod(Lon_vec, 360);
Lat_vec = max(-86, min(90, Lat_vec));

fprintf('   正在计算绝对盐度 (SA) 和守恒温度 (CT)...\n');

% 分批处理
n_valid = sum(valid_grid(:));
batch_size = 100000;
n_batches = ceil(n_valid / batch_size);

SA_vec = zeros(n_valid, 1);
CT_vec = zeros(n_valid, 1);

% 边界约束
Lon_vec = max(0.5, min(359.5, Lon_vec));
P_vec = max(0.1, min(6131, P_vec));

for batch = 1:n_batches
    idx_start = (batch-1) * batch_size + 1;
    idx_end = min(batch * batch_size, n_valid);
    
    SA_vec(idx_start:idx_end) = gsw_SA_from_SP(...
        S_vec(idx_start:idx_end), ...
        P_vec(idx_start:idx_end), ...
        Lon_vec(idx_start:idx_end), ...
        Lat_vec(idx_start:idx_end));
    
    CT_vec(idx_start:idx_end) = gsw_CT_from_pt(...
        SA_vec(idx_start:idx_end), ...
        T_vec(idx_start:idx_end));
    
    if mod(batch, 10) == 0 || batch == n_batches
        fprintf('   GSW 已处理 %d/%d 批次 (%.1f%%)\n', batch, n_batches, 100*batch/n_batches);
    end
end

fprintf('   正在调用 TEOS-10 引擎计算高阶导数...\n');
try
    Derivs_Raw = Engine.calculate_all_orders(double(SA_vec), double(CT_vec), double(P_vec), MaxOrder);
catch ME
    fprintf('   ❌ 引擎计算失败: %s\n', ME.message);
    fprintf('   ⚠️ 仅计算一阶项，高阶项置零\n');
    rho_0_local = gsw_rho(SA_vec, CT_vec, P_vec);
    alpha_local = gsw_alpha(SA_vec, CT_vec, P_vec);
    Derivs_Raw = struct();
    Derivs_Raw.d1_T = -rho_0_local .* alpha_local;
    for ord = 2:MaxOrder
        Derivs_Raw.(sprintf('d%d_T', ord)) = zeros(size(rho_0_local));
    end
end

% 将导数还原回 3D 网格
Derivs = cell(MaxOrder, 1);
for n = 1:MaxOrder
    tmp = nan(size(P_3D));
    tmp(valid_grid) = Derivs_Raw.(sprintf('d%d_T', n));
    Derivs{n} = tmp;
end

% 存储参考态 CT (用于后续异常计算)
CT_Ref = nan(size(P_3D));
CT_Ref(valid_grid) = CT_vec;

%% ========================================================================
%% 6. Phase 3: 计算 TSLA 时间序列
%% ========================================================================
fprintf('>> [5/6] Phase 3 - 计算海平面异常 (TSLA) 各阶贡献...\n');

[Nx, Ny, ~] = size(P_3D);
TotalSteps = length(Years) * length(Months);

% 预分配大矩阵
TSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder);
Time_Axis = zeros(TotalSteps, 1);

% 预计算阶乘倒数
fact_inv = zeros(1, MaxOrder);
for k = 1:MaxOrder, fact_inv(k) = 1/factorial(k); end

idx = 0;
warning('off', 'all');

for y = Years
    for m = Months
        idx = idx + 1;
        Time_Axis(idx) = y + (m-0.5)/12; % 时间中心点
        
        PatT = sprintf('*Temp*year_%d_month_%02d.nc', y, m);
        PatS = sprintf('*Salinity*year_%d_month_%02d.nc', y, m);
        
        % 文件搜索
        dT = dir(fullfile(TempDir, PatT));
        dS = dir(fullfile(SaltDir, PatS));
        
        if isempty(dT) || isempty(dS)
            TSLA_AllOrders(:,:,idx,:) = NaN;
            continue;
        end
        
        try
            % 读取瞬时场
            T_inst = double(ncread(fullfile(dT(1).folder, dT(1).name), TVar));
            S_inst = double(ncread(fullfile(dS(1).folder, dS(1).name), SVar));
            
            % [维度修正]
            if size(T_inst, 1) < size(T_inst, 2) && size(T_inst, 1) < size(T_inst, 3)
                T_inst = permute(T_inst, [2, 3, 1]);
                S_inst = permute(S_inst, [2, 3, 1]);
            end
            
            % 深度截取
            if size(T_inst, 3) > Nz
                T_inst = T_inst(:,:,1:Nz);
                S_inst = S_inst(:,:,1:Nz);
            end
            
            % 异常过滤
            T_inst(T_inst > 100 | T_inst < -10) = NaN;
            S_inst(S_inst > 100 | S_inst < 0) = NaN;
            
            % [坐标修正] 确保经纬度在 GSW 合法范围内
            LON_3D_fix = mod(LON_3D, 360);
            LAT_3D_fix = max(-90, min(90, LAT_3D));
            
            % 必须放弃近似，使用严格转换！(与 Exact 保持绝对一致)
            SA_inst = gsw_SA_from_SP(S_inst, P_3D, LON_3D_fix, LAT_3D_fix);
            CT_inst = gsw_CT_from_pt(SA_inst, T_inst);
            
            % 计算温度异常 delta_theta
            dTheta = CT_inst - CT_Ref;
            
            % 泰勒展开项计算
            for n = 1:MaxOrder
                % 公式: Term_n = (1/n!) * (Dn_rho) * (delta_theta)^n
                term = fact_inv(n) .* Derivs{n} .* (dTheta .^ n);
                
                % 垂直积分: TSLA = - (1/rho0) * ∫(rho_prime) dz * 1000 (mm)
                integral_val = nansum(term .* dz_3D, 3);
                TSLA_AllOrders(:,:,idx,n) = -(integral_val / rho0) * 1000;
            end
        catch
            TSLA_AllOrders(:,:,idx,:) = NaN;
        end
    end
    fprintf('   ... 年份 %d TSLA 计算完成\n', y);
end
warning('on', 'all');

%% ========================================================================
%% 7. 结果保存
%% ========================================================================
fprintf('>> [6/6] 保存计算结果...\n');

lat = Lat;
lon = Lon;
OutFile = fullfile(OutputDir, 'IAP_TSLA_Terms_1to8_Average.mat');
save(OutFile, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

fprintf('\n>>> ✅ IAP 数据处理完全结束! <<<\n');
fprintf('   输出文件: %s\n', OutFile);
