%% Ishii 数据处理脚本
clear; clc;
addpath('D:\work');

%% 1. 基础参数与路径配置
DataRoot = 'D:\work\Ishii_05_24';
TempDir = fullfile(DataRoot, 'Temperature');
SaltDir = fullfile(DataRoot, 'Salinity');
OutputDir = 'D:\work\Ishii_TSLA_Terms';

Years = 2005:2024;
MaxOrder = 8;
MaxDepth = 2000;
rho0 = single(1035.0);

if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 初始化高阶核心引擎
fprintf('>> [1/6] 初始化 TEOS-10 升级版热力学引擎...\n');
try
    Engine = TEOS10_General_Engine(MaxOrder);
catch
    error('❌ 无法加载 TEOS10_General_Engine。');
end

%% 3. 提取网格信息（使用2005年文件作为样本）
fprintf('>> [2/6] 提取网格信息...\n');
sample_file = fullfile(TempDir, 'temp.2005.nc');
if ~exist(sample_file, 'file'), error('❌ 找不到样本文件'); end

% 读取坐标
Lon = single(ncread(sample_file, 'longitude'));
Lat = single(ncread(sample_file, 'latitude'));
Level = single(ncread(sample_file, 'level'));  % 注意：这是 level，不是 depth

% 统一 lon 坐标
Lon = mod(Lon, 360);

% 显示坐标信息
fprintf('  经度: %d 点 (%.1f 到 %.1f)\n', length(Lon), min(Lon), max(Lon));
fprintf('  纬度: %d 点 (%.1f 到 %.1f)\n', length(Lat), min(Lat), max(Lat));
fprintf('  层次: %d 层\n', length(Level));

% 注意：Ishii 的 level 可能是压力或深度，需要确认
% 这里假设 level 是深度（米），但需要验证
fprintf('  层次值: %s\n', mat2str(Level));

%% 4. 确定有效的深度层
% 由于 Ishii 只有 28 层，我们可能需要使用所有层
DepthIdx = 1:length(Level);  % 使用所有层
CalcDepth = Level(DepthIdx);

% 检查是否有深度超过 MaxDepth 的层
deep_layers = find(CalcDepth > MaxDepth);
if ~isempty(deep_layers)
    fprintf('  警告: %d 层深度超过 %d m，将被排除\n', length(deep_layers), MaxDepth);
    DepthIdx = find(CalcDepth <= MaxDepth);
    CalcDepth = Level(DepthIdx);
end

if length(CalcDepth) < 2
    error('❌ 有效深度层数小于2，无法计算积分权重！');
end

%% 5. 生成网格
fprintf('>> [3/6] 生成网格...\n');
% Ishii 数据维度: [longitude, latitude, level, time]
% 目标处理维度: [latitude, depth, longitude]
[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Ny, Nz, Nx] = size(P_3D);

fprintf('  网格维度: %d (lat) × %d (depth) × %d (lon)\n', Ny, Nz, Nx);

% GSW 安全坐标
LON_GSW = max(0.5, min(359.5, LON_3D));
LAT_GSW = max(-86, min(90, LAT_3D));
P_GSW   = max(0.1, min(6131, P_3D));

% 计算积分权重 dz
num_depth = length(CalcDepth);
dz = zeros(num_depth, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:num_depth-1
    dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1));
end
dz(num_depth) = CalcDepth(num_depth) - CalcDepth(num_depth-1);
dz_3D = reshape(dz, 1, num_depth, 1);

%% 6. 计算多年平均态
fprintf('>> [4