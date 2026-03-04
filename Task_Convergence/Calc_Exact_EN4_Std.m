%% Calc_Formula11_Exact_Std.m
% =========================================================================
% 功能：计算【标准参考态】下的精确热容海平面 (Exact TSLA relative to Std Ref State)
% 
% 物理定义：
%   TSLA_Exact_Std = - Integral [ (rho(S_inst, T_inst, p) - rho(S_std, T_std, p)) / rho0 ] dz
% 
% 核心特征：
%   1. 参考态：S_std = 35 psu, T_std = 0 degC (常数)
%   2. 真实态：S_inst, T_inst (实时数据)
%   3. 无截断误差：直接调用 gsw_rho 函数。
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; 
rho0 = 1035.0;       
MaxDepth = 2000;     

% 标准参考态定义
Std_S_Val = 35.0;
Std_T_Val = 0.0; % degC

%% 2. 初始化网格
fprintf('>> [1/4] 初始化网格...\n');
% 自动找一个文件读取网格
SampleFile = '';
for y = Years
    List = dir(fullfile(DataDir, sprintf('*%d*',y), '*.nc'));
    if ~isempty(List)
        SampleFile = fullfile(List(1).folder, List(1).name); break; 
    end
end
if isempty(SampleFile), error('❌ 找不到数据文件'); end

Lon = ncread(SampleFile, 'lon');
Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);

% 生成 3D 坐标矩阵
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D); % 压力矩阵
[Nx, Ny, Nz] = size(P_3D);

% 计算垂直积分权重 dz
dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []); %用于广播乘法

%% 3. 预计算标准参考态密度 (rho_Std)
fprintf('>> [2/4] 预计算标准态密度 (rho_Std: S=35, T=0)...\n');

% 构造标准态 3D 矩阵
Std_S = ones(Nx, Ny, Nz) * Std_S_Val;
Std_T = ones(Nx, Ny, Nz) * Std_T_Val;


SA_Std = gsw_SA_from_SP(Std_S, P_3D, LON_3D, LAT_3D);
CT_Std = gsw_CT_from_pt(SA_Std, Std_T);
rho_Std = gsw_rho(SA_Std, CT_Std, P_3D);

%% 4. 遍历计算
fprintf('>> [3/4] 开始计算精确热容海平面 (Exact Std)...\n');

TotalSteps = length(Years) * 12;
TSLA_Exact_Std = zeros(Nx, Ny, TotalSteps, 'single'); % 结果矩阵
time_vec = zeros(TotalSteps, 1);

idx = 0;
for y = Years
    fprintf('   正在处理: %d ... ', y);
    tic;
    
    Folder = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y));
    
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        
        % 找文件
        Pattern = sprintf('*%d%02d*.nc', y, m);
        FList = dir(fullfile(Folder, Pattern));
        if isempty(FList)
            TSLA_Exact_Std(:,:,idx) = NaN; continue; 
        end
        FullPath = fullfile(FList(1).folder, FList(1).name);
        
        try
            % 读取 T, S
            T_Raw = double(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
            S_Raw = double(ncread(FullPath, 'salinity',    [1 1 1 1], [Inf Inf Nz 1]));
            
            % Kelvin 转换
            if nanmean(T_Raw(:)) > 100, T_Raw = T_Raw - 273.15; end
            
            % === 核心计算 ===
            % 1. 坐标转换
            SA_Inst = gsw_SA_from_SP(S_Raw, P_3D, LON_3D, LAT_3D);
            CT_Inst = gsw_CT_from_pt(SA_Inst, T_Raw);
            
            % 2. 计算真实密度
            rho_Inst = gsw_rho(SA_Inst, CT_Inst, P_3D);
            
            % 3. 密度差 (相对于标准态)
            drho = rho_Inst - rho_Std;
            
            % 4. 积分得到 TSLA (单位 mm)
            % steric = - integral( drho / rho0 ) dz
            val_mm = -(nansum(drho .* dz_perm, 3) / rho0) * 1000;
            
            % 屏蔽无效点
            Mask = isnan(T_Raw(:,:,1));
            val_mm(Mask) = NaN;
            
            TSLA_Exact_Std(:,:,idx) = single(val_mm);
            
        catch
            TSLA_Exact_Std(:,:,idx) = NaN;
        end
    end
    fprintf('耗时 %.1f 秒\n', toc);
end

%% 5. 保存结果
fprintf('>> [4/4] 保存文件...\n');
SaveName = fullfile(OutputDir, 'EN4_Formula11_Exact_Std.mat');
lon = Lon; lat = Lat;
save(SaveName, 'TSLA_Exact_Std', 'lon', 'lat', 'time_vec', '-v7.3');

fprintf('>> ✅ 计算完成！\n');
fprintf('   文件已生成: %s\n', SaveName);
