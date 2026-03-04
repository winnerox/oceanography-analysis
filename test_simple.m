%% 简单测试脚本
clear; clc;
addpath('D:\work');

%% 检查文件是否存在
DataRoot = 'D:\work\IAP_05_24';
TempDir = fullfile(DataRoot, 'TEMP');
SampleFile = fullfile(TempDir, 'IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc');

if ~exist(SampleFile, 'file')
    error('❌ 找不到样本文件: %s', SampleFile);
end

%% 读取网格信息
fprintf('读取网格信息...\n');
Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, 'depth_std'));

fprintf('Lon 维度: %d\n', length(Lon));
fprintf('Lat 维度: %d\n', length(Lat));
fprintf('Depth 维度: %d\n', length(Depth));
fprintf('Depth 值: %s\n', mat2str(Depth(1:min(10, length(Depth)))));

%% 检查深度
MaxDepth = 2000;
DepthIdx = find(Depth <= MaxDepth);
fprintf('深度 <= %d 的索引数量: %d\n', MaxDepth, length(DepthIdx));
if isempty(DepthIdx)
    error('❌ 没有找到深度 <= %d 的点！', MaxDepth);
end

CalcDepth = Depth(DepthIdx);
fprintf('CalcDepth 长度: %d\n', length(CalcDepth));
fprintf('CalcDepth 值: %s\n', mat2str(CalcDepth));

if length(CalcDepth) < 2
    error('❌ CalcDepth 长度小于2，无法计算积分权重！');
end

%% 测试积分权重计算
num_depth = length(CalcDepth);
dz = zeros(num_depth, 1, 'single');
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:num_depth-1
    dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1));
end
dz(num_depth) = CalcDepth(num_depth) - CalcDepth(num_depth-1);

fprintf('积分权重计算成功！\n');
fprintf('dz 值: %s\n', mat2str(dz));

%% 测试网格生成
[LAT_3D, LON_3D, DEPTH_3D] = ndgrid(Lat, Lon, CalcDepth);
fprintf('网格维度: %d x %d x %d\n', size(LAT_3D, 1), size(LAT_3D, 2), size(LAT_3D, 3));

fprintf('✅ 所有测试通过！\n');