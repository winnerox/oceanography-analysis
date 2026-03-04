%%Calc_Formula11_Exact.m
% =========================================================================
% 功能：计算精确热容海平面 (Exact Thermosteric Sea Level) - Formula 11
% 
% 物理定义：
%   TSLA_Exact = Integral [ 1 - rho(S, T, p) / rho(S, 0, p) ] dz
%   或者等价于密度差积分：
%   TSLA_Exact ≈ - Integral [ (rho(S,T) - rho(S,0)) / rho0 ] dz
% 
% 核心特征：
%   1. 控制变量：盐度始终使用“实时真实盐度” (S_inst)。
%   2. 变温过程：温度从 0°C (参考态) 变到 T_inst (真实态)。
%   3. 无截断误差：直接调用 gsw_rho 函数，不使用泰勒展开。
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; 
rho0 = 1035.0;       
MaxDepth = 2000;     
Std_Temp_Val = 0.0;  % 参考态位温 = 0度

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

%% 3. 遍历计算
fprintf('>> [2/4] 开始计算精确热容海平面 (No Taylor Expansion)...\n');

TotalSteps = length(Years) * 12;
TSLA_Thermo_Exact = zeros(Nx, Ny, TotalSteps, 'single'); % 结果矩阵
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
            TSLA_Thermo_Exact(:,:,idx) = NaN; continue; 
        end
        FullPath = fullfile(FList(1).folder, FList(1).name);
        
        try
            % 读取 T, S
            T_Raw = double(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
            S_Raw = double(ncread(FullPath, 'salinity',    [1 1 1 1], [Inf Inf Nz 1]));
            
            % Kelvin 转换
            if nanmean(T_Raw(:)) > 100, T_Raw = T_Raw - 273.15; end
            
            % === 核心计算 (精确差值法) ===
            % 1. 坐标转换
            SA = gsw_SA_from_SP(S_Raw, P_3D, LON_3D, LAT_3D);
            CT_Inst = gsw_CT_from_pt(SA, T_Raw);
            
            % 2. 构造参考态温度 (T=0, S=Real)
            % 注意：gsw_CT_from_pt 要求第二个参数也是矩阵
            CT_Ref = gsw_CT_from_pt(SA, Std_Temp_Val * ones(size(SA)));
            
            % 3. 直接调用 GSW 算密度 (这是最准的，包含所有非线性)
            rho_Inst = gsw_rho(SA, CT_Inst, P_3D); % 真实状态
            rho_Ref  = gsw_rho(SA, CT_Ref,  P_3D); % 参考状态 (仅T变为0)
            
            % 4. 密度差
            drho = rho_Inst - rho_Ref;
            
            % 5. 积分得到 TSLA (单位 mm)
            % 公式: steric = - integral( drho / rho0 ) dz
            val_mm = -(nansum(drho .* dz_perm, 3) / rho0) * 1000;
            
            % 屏蔽无效点
            Mask = isnan(T_Raw(:,:,1));
            val_mm(Mask) = NaN;
            
            TSLA_Thermo_Exact(:,:,idx) = single(val_mm);
            
        catch
            TSLA_Thermo_Exact(:,:,idx) = NaN;
        end
    end
    fprintf('耗时 %.1f 秒\n', toc);
end

%% 4. 保存结果
fprintf('>> [3/4] 保存文件...\n');
SaveName = fullfile(OutputDir, 'EN4_Formula11_Exact_Thermal.mat');
lon = Lon; lat = Lat;
save(SaveName, 'TSLA_Thermo_Exact', 'lon', 'lat', 'time_vec', '-v7.3');

fprintf('>> ✅ 计算完成！\n');
fprintf('   文件已生成: %s\n', SaveName);
fprintf('   现在你可以去运行 Step9 进行最终对比了。\n');