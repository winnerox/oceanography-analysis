%% EN4_TSLA_Average.m
% =========================================================================
% EN4 (Met Office) 海洋数据 TSLA 计算脚本 [平均态版本]
% =========================================================================
% 功能说明：
%   基于 Met Office Hadley Centre EN4.2.2 数据集计算 TSLA。
%   计算相对于 2005-2024 年平均态的热比容海平面异常。
%
% 核心特征：
%   1. 数据格式：单月单文件 (EN.4.2.2.f.analysis.c13.YYYYMM.nc)。
%   2. 物理量：温度为开尔文 (Kelvin)，需减去 273.15 转换为摄氏度。
%   3. 深度校正：EN4 深度层不均匀，需正确应用权重。
%
% 输入数据：
%   - EN4_analyses_c13_last20years\EN.4.2.2.analyses.c13.YYYY\EN.4.2.2.f.analysis.c13.YYYYMM.nc
%
% 输出数据：
%   - EN4_TSLA_Terms_1to8_Average.mat
% =========================================================================
clear; clc;

%% ========================================================================
%% 1. 基础参数与路径配置
%% ========================================================================
DataRoot = 'D:\work\EN4_analyses_c13_last20years';
OutputDir = 'D:\work\MAT_Data';

Years = 2005:2024;
Months = 1:12;

% 计算参数
MaxOrder = 8;
MaxDepth = 2000;
rho0 = 1035.0;

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

% 3.1 搜索样本文件
SampleFile = '';
for y = Years
    YearDir = fullfile(DataRoot, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = 1:12
        Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        FileList = dir(fullfile(YearDir, Pattern));
        if ~isempty(FileList)
            SampleFile = fullfile(FileList(1).folder, FileList(1).name);
            break;
        end
    end
    if ~isempty(SampleFile), break; end
end

if isempty(SampleFile), error('❌ 未找到 EN4 样本文件'); end
fprintf('   样本文件: %s\n', SampleFile);

% 3.2 读取元数据
warning('off', 'all');
Info = ncinfo(SampleFile);
warning('on', 'all');

VarNames = {Info.Variables.Name};
fprintf('   检测变量: %s\n', strjoin(VarNames, ', '));

% 变量映射
if any(strcmpi('temperature', VarNames)), TVar = 'temperature';
else, error('❌ 未找到 temperature 变量'); end

if any(strcmpi('salinity', VarNames)), SVar = 'salinity';
else, error('❌ 未找到 salinity 变量'); end

LonVar = 'lon'; LatVar = 'lat'; DepthVar = 'depth'; % EN4 标准命名

% 3.3 读取坐标
Lon = double(ncread(SampleFile, LonVar));
Lat = double(ncread(SampleFile, LatVar));
Depth = double(ncread(SampleFile, DepthVar));

DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
Nz = length(CalcDepth);
fprintf('   计算网格: %d x %d x %d\n', length(Lon), length(Lat), Nz);

% 3.4 构建 3D 网格
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);

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
%% 4. Phase 1: 计算多年平均态
%% ========================================================================
fprintf('>> [3/6] Phase 1 - 计算 20 年平均气候态...\n');

Sum_T = zeros(size(P_3D));
Sum_S = zeros(size(P_3D));
Count = zeros(size(P_3D), 'uint16');

warning('off', 'all');
for y = Years
    YearDir = fullfile(DataRoot, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = Months
        Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        d = dir(fullfile(YearDir, Pattern));
        if isempty(d), continue; end
        
        File = fullfile(d(1).folder, d(1).name);
        try
            % 维度: [Lon, Lat, Depth, Time=1]
            T_raw = double(ncread(File, TVar));
            S_raw = double(ncread(File, SVar));
            
            % 提取第一维度作为 Lon (EN4 通常正确，无需 permute)
            if size(T_raw, 3) > Nz
                T_raw = T_raw(:,:,1:Nz);
                S_raw = S_raw(:,:,1:Nz);
            end
            
            % 【关键】开尔文转摄氏度
            T_raw = T_raw - 273.15;
            
            % 异常值过滤
            T_raw(T_raw > 100 | T_raw < -10) = NaN;
            S_raw(S_raw > 100 | S_raw < 0) = NaN;
            
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

Mean_T = Sum_T ./ double(Count);
Mean_S = Sum_S ./ double(Count);
Mean_T(Count == 0) = NaN;
Mean_S(Count == 0) = NaN;

clear Sum_T Sum_S Count;

%% ========================================================================
%% 5. Phase 2: 计算背景场热力学导数
%% ========================================================================
fprintf('>> [4/6] Phase 2 - 计算背景场高阶导数...\n');

valid_grid = ~isnan(Mean_T) & ~isnan(Mean_S);
S_vec = double(Mean_S(valid_grid));
T_vec = double(Mean_T(valid_grid));
P_vec = double(P_3D(valid_grid));
Lon_vec = double(LON_3D(valid_grid));
Lat_vec = double(LAT_3D(valid_grid));

Lon_vec = mod(Lon_vec, 360);
Lat_vec = max(-90, min(90, Lat_vec));

try
    SA_vec = gsw_SA_from_SP(S_vec, P_vec, Lon_vec, Lat_vec);
    CT_vec = gsw_CT_from_pt(SA_vec, T_vec);
catch ME
    fprintf('   ⚠️ GSW 近似模式: %s\n', ME.message);
    SA_vec = S_vec * 1.0047;
    CT_vec = T_vec;
end

try
    Derivs_Raw = Engine.calculate_all_orders(double(SA_vec), double(CT_vec), double(P_vec), MaxOrder);
catch ME
    fprintf('   ❌ 引擎降级: %s\n', ME.message);
    rho_0_l = gsw_rho(SA_vec, CT_vec, P_vec);
    alpha_l = gsw_alpha(SA_vec, CT_vec, P_vec);
    Derivs_Raw = struct('d1_T', -rho_0_l .* alpha_l);
    for ord=2:MaxOrder, Derivs_Raw.(sprintf('d%d_T', ord)) = zeros(size(rho_0_l)); end
end

Derivs = cell(MaxOrder, 1);
for n = 1:MaxOrder
    tmp = nan(size(P_3D));
    tmp(valid_grid) = Derivs_Raw.(sprintf('d%d_T', n));
    Derivs{n} = tmp;
end

CT_Ref = nan(size(P_3D));
if exist('CT_vec', 'var'), CT_Ref(valid_grid) = CT_vec;
else, CT_Ref(valid_grid) = T_vec; end

%% ========================================================================
%% 6. Phase 3: 计算 TSLA 时间序列
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
    YearDir = fullfile(DataRoot, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = Months
        idx = idx + 1;
        Time_Axis(idx) = y + (m-0.5)/12;
        
        Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        d = dir(fullfile(YearDir, Pattern));
        if isempty(d)
            TSLA_AllOrders(:,:,idx,:) = NaN;
            continue;
        end
        
        File = fullfile(d(1).folder, d(1).name);
        try
            T_inst = double(ncread(File, TVar));
            S_inst = double(ncread(File, SVar));
            
            if size(T_inst, 3) > Nz
                T_inst = T_inst(:,:,1:Nz);
                S_inst = S_inst(:,:,1:Nz);
            end
            
            % 开尔文转摄氏度
            T_inst = T_inst - 273.15;
            
            T_inst(T_inst > 100 | T_inst < -10) = NaN;
            S_inst(S_inst > 100 | S_inst < 0) = NaN;
            
            % 正确转换：势温 -> 保守温度
            SA_inst = gsw_SA_from_SP(S_inst, P_3D, LON_3D, LAT_3D);
            CT_inst = gsw_CT_from_pt(SA_inst, T_inst);
            
            dTheta = CT_inst - CT_Ref;
            
            for n = 1:MaxOrder
                term = fact_inv(n) .* Derivs{n} .* (dTheta .^ n);
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
fprintf('>> [6/6] 保存结果...\n');

lat = Lat;
lon = Lon;
OutFile = fullfile(OutputDir, 'EN4_TSLA_Terms_1to8_Average.mat');
save(OutFile, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

fprintf('\n>>> ✅ EN4 数据处理完全结束! <<<\n');
fprintf('   输出文件: %s\n', OutFile);