%% IAP_TSLA_StdRef.m
% =========================================================================
% IAP (CAS) 海洋数据 TSLA 计算脚本 [标准态版本]
% =========================================================================
% 核心对齐：完全采用与 Calc_S_Terms_IAP_Std.m 一致的稳健架构
% =========================================================================
clear; clc;
addpath('D:\work');

%% 1. 基础参数与路径配置
DataRoot = 'D:\work\IAP_05_24';
TempDir = fullfile(DataRoot, 'TEMP');
SaltDir = fullfile(DataRoot, 'SALT');
OutputDir = 'D:\work\MAT_Data';

Years = 2005:2024;
Months = 1:12;

MaxOrder = 8;
MaxDepth = 2000;
rho0 = single(1035.0);

Std_Salt_Val = single(35.0);  
Std_Temp_Val = single(0.0);   

if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 初始化核心引擎 (使用高阶)
fprintf('>> [1/6] 初始化 TEOS-10 高阶热力学引擎...\n');
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
try
    Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);
catch
    error('❌ 无法加载 TEOS10_HighOrder_Engine。');
end

%% 3. 网格与坐标 (对齐 Calc_S_Terms_IAP_Std)
fprintf('>> [2/6] 正在侦测网格与变量定义...\n');
List = dir(fullfile(TempDir, '*.nc'));
if isempty(List), error('❌ 未找到样本文件'); end
SampleFile = fullfile(List(1).folder, List(1).name);

vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth_std')), vDepth='depth_std';
elseif any(strcmp(vnames, 'depth')), vDepth='depth'; 
else, vDepth='depth_std'; end

if any(strcmpi('temp', vnames)), TVar = 'temp'; else, TVar = 'temperature'; end

SampleS_File = strrep(SampleFile, 'TEMP', 'SALT');
SampleS_File = strrep(SampleS_File, 'Temp', 'Salinity');
SampleS_File = strrep(SampleS_File, 'IAPv4', 'IAPv2');
vi_s = ncinfo(SampleS_File); vnames_s = {vi_s.Variables.Name};
if any(strcmpi('salinity', vnames_s)), SVar = 'salinity'; else, SVar = 'sal'; end

% 采用原生读取，不做任何坐标破坏！
Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, vDepth)); 
DepthIdx = find(Depth <= MaxDepth); 
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);

% 边界约束
% 【关键】GSW SAAR 参考数据纬度范围是 [-86, 90]，超出范围会导致索引错误
LON_3D = max(0.5, min(359.5, LON_3D));
LAT_3D = max(-86, min(90, LAT_3D));  % 限制在 GSW 有效范围内

% 压力场计算保持原生
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
P_3D = max(0.1, min(6131, P_3D));    % GSW 压力范围 [0, 6131]
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_3D = reshape(dz, 1, 1, []);

%% 4. Phase 1: 构建标准态参考场
fprintf('>> [3/6] Phase 1 - 构建标准态参考场...\n');

% 与 S_Terms 完全一致的强制类型转换，绕开 GSW 查表 Bug
try
    SA_Ref_3D = single(gsw_SA_from_SP(double(Std_Salt_Val * ones(Nx, Ny, Nz, 'single')), double(P_3D), double(LON_3D), double(LAT_3D)));
    CT_Ref_3D = single(gsw_CT_from_pt(SA_Ref_3D, double(Std_Temp_Val * ones(Nx, Ny, Nz, 'single'))));
catch ME
    error('❌ GSW 工具箱报错: %s', ME.message);
end

%% 5. Phase 2: 计算高阶纯 T 导数
fprintf('>> [4/6] Phase 2 - 计算标准态下的高阶导数...\n');

Derivs_T = cell(MaxOrder, 1);
try
    for n = 1:MaxOrder
        % 计算纯温度项导数 (n_T = n, n_S = 0)
        D_val = Engine.calculate_mixed(double(SA_Ref_3D), double(CT_Ref_3D), double(P_3D), n, 0);
        Derivs_T{n} = single(D_val);
    end
catch ME
    error('❌ 高阶热力学引擎计算失败: %s', ME.message);
end

InvFact = zeros(21,1,'single'); 
for i=0:20, InvFact(i+1) = single(1/factorial(i)); end

%% 6. Phase 3: 计算 TSLA 时间序列
fprintf('>> [5/6] Phase 3 - 计算海平面异常 (TSLA)...\n');
TotalSteps = length(Years) * length(Months);
TSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
Time_Axis = zeros(TotalSteps, 1);
idx = 0; warning('off', 'all');

for y = Years
    for m = Months
        idx = idx + 1;
        Time_Axis(idx) = y + (m-0.5)/12;
        
        PatT = sprintf('*Temp*year_%d_month_%02d.nc', y, m);
        PatS = sprintf('*Salinity*year_%d_month_%02d.nc', y, m);
        dT = dir(fullfile(TempDir, PatT)); dS = dir(fullfile(SaltDir, PatS));
        
        if isempty(dT) || isempty(dS)
            TSLA_AllOrders(:,:,idx,:) = NaN; continue;
        end
        
        try
            FileT = fullfile(dT(1).folder, dT(1).name);
            FileS = fullfile(dS(1).folder, dS(1).name);
            
            % 读取并转维
            T_inst = single(permute(ncread(FileT, TVar, [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
            S_inst = single(permute(ncread(FileS, SVar, [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
            
            % =================================================================
            % 🚀 性能加速核心：提取一维真实海洋点，绝不给陆地算热力学！
            % =================================================================
            % 找到有效水体的布尔掩膜
            MaskValid = ~(isnan(T_inst) | isnan(S_inst) | (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0));
            idx_valid = find(MaskValid); % 提取一维索引
            
            % 仅把真实海水的 1D 向量（约占全体积的 30%）转为 double 投喂给 GSW
            SA_1D = single(gsw_SA_from_SP(double(S_inst(idx_valid)), double(P_3D(idx_valid)), double(LON_3D(idx_valid)), double(LAT_3D(idx_valid))));
            CT_1D = single(gsw_CT_from_pt(SA_1D, double(T_inst(idx_valid))));
            
            % 计算异常并映射回 3D 矩阵
            dTheta = zeros(Nx, Ny, Nz, 'single');
            dTheta(idx_valid) = CT_1D - CT_Ref_3D(idx_valid);
            
            % 🚀 性能加速 2：用累乘 (.*) 替代极其耗时的幂运算 (.^ n)
            dTheta_n = dTheta; % 初始化 1 阶 (即 dTheta^1)
            
            for n = 1:MaxOrder
                C = InvFact(n+1);
                Deriv = Derivs_T{n};
                
                % 此时 dTheta_n 已经是当前的 dTheta^n 了
                rho_t = C .* Deriv .* dTheta_n;
                rho_t(~MaskValid) = 0; % 确保陆地不参与积分
                
                val = -single(sum(rho_t .* dz_3D, 3)) / rho0 * single(1000);
                
                % 还原 NaN 掩膜
                valid_depths = single(sum(MaskValid, 3));
                val(valid_depths < 1) = NaN;
                
                TSLA_AllOrders(:,:,idx,n) = val;
                
                % 为下一阶准备：自身再乘一次 dTheta (比 dTheta .^ (n+1) 快得多)
                if n < MaxOrder
                    dTheta_n = dTheta_n .* dTheta;
                end
            end
        catch ME
            fprintf('   ⚠️ 计算出错 %d-%02d: %s\n', y, m, ME.message);
            TSLA_AllOrders(:,:,idx,:) = NaN;
        end
    end
    fprintf('   ... 年份 %d 计算完成\n', y);
end
warning('on', 'all');

%% 7. 结果保存
fprintf('>> [6/6] 保存结果...\n');
lat = Lat; lon = Lon;
OutFile = fullfile(OutputDir, 'IAP_TSLA_Terms_1to8_StdRef.mat');
save(OutFile, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');
fprintf('\n>>> ✅ IAP 标准态处理完全结束! <<<\n');