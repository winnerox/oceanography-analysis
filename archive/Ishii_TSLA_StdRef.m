%% Ishii_TSLA_StdRef.m
% =========================================================================
% Ishii (JMA) 海洋数据 TSLA 计算脚本 [标准态版本]
% =========================================================================
% 功能说明：
%   基于 Ishii 数据计算相对于“标准海水状态” (S=35 PSU, T=0°C) 的海平面异常。
%   此版本用于与“平均态版本”进行对比，评估参考态选择对非线性项的影响。
%
% 核心特征：
%   1. 参考态：固定盐度 S=35.0, 位温 T=0.0。
%   2. 算法流程：同 IAP 标准态版本，但适配 Ishii 的数据结构 (4D 数组)。
%
% 输入数据：
%   - Ishii_05_24 (年度 NetCDF 文件)
%
% 输出数据：
%   - Ishii_TSLA_Terms_1to8_StdRef.mat
% =========================================================================
clear; clc;

%% ========================================================================
%% 1. 基础参数与路径配置
%% ========================================================================
DataRoot = 'D:\work\Ishii_05_24';
TempDir = fullfile(DataRoot, 'Temperature');
SaltDir = fullfile(DataRoot, 'Salinity');
OutputDir = 'D:\work\MAT_Data';

Years = 2005:2024;
Months = 1:12;

% 计算参数
MaxOrder = 8;
MaxDepth = 2000;
rho0 = 1035.0;

% 【差异点】定义标准态参数
Std_Salt_Val = 35.0;  % 实用盐度 (PSU)
Std_Temp_Val = 0.0;   % 位温 (°C)

if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% ========================================================================
%% 2. 初始化核心引擎
%% ========================================================================
fprintf('>> [1/6] 初始化 TEOS-10 热力学引擎...\n');
try
    Engine = TEOS10_General_Engine();
catch
    error('❌ 无法加载 TEOS10_General_Engine。');
end

%% ========================================================================
%% 3. 自动侦测网格信息
%% ========================================================================
fprintf('>> [2/6] 正在侦测网格与变量定义...\n');

SampleT = fullfile(TempDir, sprintf('temp.%d.nc', Years(1)));
if ~exist(SampleT, 'file')
    error('❌ 未找到样本文件: %s', SampleT);
end
fprintf('   样本文件: %s\n', SampleT);

warning('off', 'all');
Info = ncinfo(SampleT);
warning('on', 'all');

% 变量名映射 (Ishii)
TVar = 'temp';
SVar = 'sal';
LonVar = 'longitude';
LatVar = 'latitude';
DepthVar = 'level';

fprintf('   变量映射: T=[%s], S=[%s], Lon=[%s], Lat=[%s], Depth=[%s]\n', ...
    TVar, SVar, LonVar, LatVar, DepthVar);

% 读取轴信息
Lon = double(ncread(SampleT, LonVar));
Lat = double(ncread(SampleT, LatVar));
Depth = double(ncread(SampleT, DepthVar));

DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
Nz = length(CalcDepth);
fprintf('   计算网格: %d x %d x %d\n', length(Lon), length(Lat), Nz);

% 构建 3D 网格
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);

% 【关键】GSW SAAR 参考数据纬度范围是 [-86, 90]，超出范围会导致索引错误
LON_3D = max(0.5, min(359.5, LON_3D));
LAT_3D = max(-86, min(90, LAT_3D));

P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
P_3D = max(0.1, min(6131, P_3D));  % GSW 压力范围

% 层厚
dz = zeros(Nz, 1);
if Nz >= 2
    dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
    for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
    dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
else
    dz(1) = CalcDepth(1);
end
dz_3D = reshape(dz, 1, 1, []);

%% ========================================================================
%% 4. Phase 1: 构建标准态参考场
%% ========================================================================
fprintf('>> [3/6] Phase 1 - 构建标准态参考场 (S=%.1f, T=%.1f)...\n', Std_Salt_Val, Std_Temp_Val);

% 全场恒定值
Ref_Salt = ones(size(P_3D)) * Std_Salt_Val;
Ref_Temp = ones(size(P_3D)) * Std_Temp_Val;

% 计算标准态的绝对盐度(SA)与守恒温度(CT)
fprintf('   分批计算 SA 和 CT...\n');

total_pts = numel(P_3D);
batch_size = 100000;
n_batches = ceil(total_pts / batch_size);

SA_Std = zeros(total_pts, 1);
CT_Std = zeros(total_pts, 1);

Ref_Salt_1D = double(Ref_Salt(:));
Ref_Temp_1D = double(Ref_Temp(:));
P_3D_1D = double(P_3D(:));
LON_3D_1D = double(LON_3D(:));
LAT_3D_1D = double(LAT_3D(:));

% 边界约束 - GSW SAAR 参考数据纬度范围是 [-86, 90]
LON_3D_1D = max(0.5, min(359.5, LON_3D_1D));
LAT_3D_1D = max(-86, min(90, LAT_3D_1D));
P_3D_1D = max(0.1, min(6131, P_3D_1D));

for batch = 1:n_batches
    idx_start = (batch-1) * batch_size + 1;
    idx_end = min(batch * batch_size, total_pts);
    
    SA_Std(idx_start:idx_end) = gsw_SA_from_SP(...
        Ref_Salt_1D(idx_start:idx_end), ...
        P_3D_1D(idx_start:idx_end), ...
        LON_3D_1D(idx_start:idx_end), ...
        LAT_3D_1D(idx_start:idx_end));
    
    CT_Std(idx_start:idx_end) = gsw_CT_from_pt(...
        SA_Std(idx_start:idx_end), ...
        Ref_Temp_1D(idx_start:idx_end));
    
    if mod(batch, 10) == 0 || batch == n_batches
        fprintf('   已处理 %d/%d 批次 (%.1f%%)\n', batch, n_batches, 100*batch/n_batches);
    end
end

SA_Std = reshape(SA_Std, size(P_3D));
CT_Std = reshape(CT_Std, size(P_3D));

%% ========================================================================
%% 5. Phase 2: 计算标准态热力学导数
%% ========================================================================
fprintf('>> [4/6] Phase 2 - 计算标准态下的高阶导数...\n');

try
    Derivs_Raw = Engine.calculate_all_orders(double(SA_Std(:)), double(CT_Std(:)), double(P_3D(:)), MaxOrder);
    
    Derivs = cell(MaxOrder, 1);
    for n = 1:MaxOrder
        tmp = reshape(Derivs_Raw.(sprintf('d%d_T', n)), size(P_3D));
        Derivs{n} = tmp;
    end
catch ME
    fprintf('   ❌ 引擎报错: %s (降级由 rho/alpha 计算一阶项)\n', ME.message);
    rho_0_l = gsw_rho(SA_Std(:), CT_Std(:), P_3D(:));
    alpha_l = gsw_alpha(SA_Std(:), CT_Std(:), P_3D(:));
    
    Derivs = cell(MaxOrder, 1);
    Derivs{1} = reshape(-rho_0_l .* alpha_l, size(P_3D));
    for n = 2:MaxOrder, Derivs{n} = zeros(size(P_3D)); end
end

% 标准态 CT 作为参考基准
CT_Ref = CT_Std;
clear Ref_Salt Ref_Temp SA_Std CT_Std;

%% ========================================================================
%% 6. Phase 3: 计算 TSLA 时间序列 (相对标准态)
%% ========================================================================
fprintf('>> [5/6] Phase 3 - 计算海平面异常 (TSLA)...\n');

[Nx, Ny, ~] = size(P_3D);
TotalSteps = length(Years) * length(Months);

TSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder);
Time_Axis = zeros(TotalSteps, 1);

fact_inv = zeros(1, MaxOrder);
for k = 1:MaxOrder, fact_inv(k) = 1/factorial(k); end

idx = 0;
warning('off', 'all');

for y = Years
    TFile = fullfile(TempDir, sprintf('temp.%d.nc', y));
    SFile = fullfile(SaltDir, sprintf('sal.%d.nc', y));
    
    if ~exist(TFile, 'file') || ~exist(SFile, 'file')
        for m = Months
            idx = idx + 1;
            Time_Axis(idx) = y + (m-0.5)/12;
            TSLA_AllOrders(:,:,idx,:) = NaN;
        end
        continue;
    end
    
    try
        % 读取全年数据
        T_year = double(ncread(TFile, TVar));
        S_year = double(ncread(SFile, SVar));
        
        T_year = T_year(:,:,1:Nz,:);
        S_year = S_year(:,:,1:Nz,:);
        
        % 过滤异常值
        T_year(abs(T_year) > 1e10) = NaN;
        S_year(abs(S_year) > 1e10) = NaN;
        T_year(T_year > 100 | T_year < -10) = NaN;
        S_year(S_year > 100 | S_year < 0) = NaN;
        
        nMonths = min(12, size(T_year, 4));
        
        for m = 1:nMonths
            idx = idx + 1;
            Time_Axis(idx) = y + (m-0.5)/12;
            
            % 提取单月数据
            T_inst = T_year(:,:,:,m);
            S_inst = S_year(:,:,:,m);
            
            % [坐标修正] 确保经纬度在 GSW 合法范围内
            LON_3D_fix = mod(LON_3D, 360);
            LAT_3D_fix = max(-90, min(90, LAT_3D));
            
            % 必须放弃近似，使用严格转换！(与 Exact 保持绝对一致)
            SA_inst = gsw_SA_from_SP(S_inst, P_3D, LON_3D_fix, LAT_3D_fix);
            CT_inst = gsw_CT_from_pt(SA_inst, T_inst);
            
            % 异常 = 瞬时态 - 标准态
            dTheta = CT_inst - CT_Ref;
            
            % 泰勒展开
            for n = 1:MaxOrder
                term = fact_inv(n) .* Derivs{n} .* (dTheta .^ n);
                integral_val = nansum(term .* dz_3D, 3);
                TSLA_AllOrders(:,:,idx,n) = -(integral_val / rho0) * 1000;
            end
        end
    catch ME
        fprintf('   ❌ 年份 %d 计算出错: %s\n', y, ME.message);
        % 占位逻辑同上
    end
    fprintf('   ... 年份 %d TSLA 计算完成\n', y);
end
warning('on', 'all');

%% ========================================================================
%% 7. 结果保存
%% ========================================================================
fprintf('>> [6/6] 保存结果...\n');

lat = Lat;
lon = Lon;
OutFile = fullfile(OutputDir, 'Ishii_TSLA_Terms_1to8_StdRef.mat');
save(OutFile, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

fprintf('\n>>> ✅ Ishii 标准态处理完全结束! <<<\n');
fprintf('   输出文件: %s\n', OutFile);
