%% EN4_TSLA_StdRef.m
% =========================================================================
% EN4 (Met Office) 海洋数据 TSLA 计算脚本 [标准态版本]
% =========================================================================
% 功能说明：
%   基于 EN4 数据计算相对于“标准海水状态” (S=35 PSU, T=0°C) 的海平面异常。
%   此版本用于与“平均态版本”进行对比。
%
% 核心特征：
%   1. 参考态：S=35.0, T=0.0。
%   2. 温度单位：需处理 Kelvin -> Celsius 转换。
%
% 输入数据：
%   - EN4_analyses_c13_last20years\EN.4.2.2.analyses.c13.YYYY\EN.4.2.2.f.analysis.c13.YYYYMM.nc
%
% 输出数据：
%   - EN4_TSLA_Terms_1to8_StdRef.mat
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

% 【标准态】
Std_Salt_Val = 35.0;
Std_Temp_Val = 0.0;

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
%% 3. 自动侦测网格信息 (EN4)
%% ========================================================================
fprintf('>> [2/6] 正在侦测网格与变量定义...\n');

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

warning('off', 'all');
Info = ncinfo(SampleFile);
warning('on', 'all');

VarNames = {Info.Variables.Name};
if any(strcmpi('temperature', VarNames)), TVar = 'temperature'; else, error('❌'); end
if any(strcmpi('salinity', VarNames)), SVar = 'salinity'; else, error('❌'); end
LonVar = 'lon'; LatVar = 'lat'; DepthVar = 'depth';

Lon = double(ncread(SampleFile, LonVar));
Lat = double(ncread(SampleFile, LatVar));
Depth = double(ncread(SampleFile, DepthVar));

DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
Nz = length(CalcDepth);
fprintf('   计算网格: %d x %d x %d\n', length(Lon), length(Lat), Nz);

[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);

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

Ref_Salt = ones(size(P_3D)) * Std_Salt_Val;
Ref_Temp = ones(size(P_3D)) * Std_Temp_Val;

LON_3D_fix = mod(LON_3D, 360);
LAT_3D_fix = max(-90, min(90, LAT_3D));

try
    SA_Std = gsw_SA_from_SP(Ref_Salt, P_3D, LON_3D_fix, LAT_3D_fix);
    CT_Std = gsw_CT_from_pt(SA_Std, Ref_Temp);
catch
    fprintf('   ⚠️ GSW 计算失败，使用简化换算...\n');
    SA_Std = Ref_Salt * 1.0047;
    CT_Std = Ref_Temp;
end

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
    fprintf('   ❌ 引擎报错: %s\n', ME.message);
    rho_0_l = gsw_rho(SA_Std(:), CT_Std(:), P_3D(:));
    alpha_l = gsw_alpha(SA_Std(:), CT_Std(:), P_3D(:));
    
    Derivs = cell(MaxOrder, 1);
    Derivs{1} = reshape(-rho_0_l .* alpha_l, size(P_3D));
    for n = 2:MaxOrder, Derivs{n} = zeros(size(P_3D)); end
end

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
            
            % Kelvin -> Celsius
            T_inst = T_inst - 273.15;
            
            T_inst(T_inst > 100 | T_inst < -10) = NaN;
            S_inst(S_inst > 100 | S_inst < 0) = NaN;
            
            % 正确转换：势温 -> 保守温度
            SA_inst = gsw_SA_from_SP(S_inst, P_3D, LON_3D, LAT_3D);
            CT_inst = gsw_CT_from_pt(SA_inst, T_inst);
            
            % 异常
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
OutFile = fullfile(OutputDir, 'EN4_TSLA_Terms_1to8_StdRef.mat');
save(OutFile, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

fprintf('\n>>> ✅ EN4 标准态处理完全结束! <<<\n');
fprintf('   输出文件: %s\n', OutFile);