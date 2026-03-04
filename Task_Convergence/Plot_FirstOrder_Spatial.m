%% Plot_FirstOrder_Spatial.m
% =========================================================================
% 功能：绘制三个数据集（EN4、IAP、Ishii）的平均态一阶项空间分布
% 输出：
%   1. 一阶项的空间分布图（3D视角）
%   2. 一阶项的空间分布图（等高线）
% =========================================================================
clear; clc;
addpath('D:\work');

% 定义数据文件路径
EN4_File = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms_1to8_Average.mat';
IAP_File = 'D:\work\IAP_mat_data\IAP_TSLA_Terms_1to8_Average.mat';
Ishii_File = 'D:\work\Ishii_mat_data\Ishii_TSLA_Terms_1to8_Average.mat';

% 加载数据
fprintf('加载数据...\n');
load(EN4_File, 'TSLA_AllOrders', 'lat', 'lon');
EN4_FirstOrder = squeeze(TSLA_AllOrders(:, :, :, 1));
EN4_lat = lat;
EN4_lon = lon;
EN4_mean = squeeze(mean(EN4_FirstOrder, 3));
EN4_mean(isnan(EN4_mean)) = 0;

load(IAP_File, 'TSLA_AllOrders', 'lat', 'lon');
IAP_FirstOrder = squeeze(TSLA_AllOrders(:, :, :, 1));
IAP_lat = lat;
IAP_lon = lon;
IAP_mean = squeeze(mean(IAP_FirstOrder, 3));
IAP_mean(isnan(IAP_mean)) = 0;

load(Ishii_File, 'TSLA_AllOrders', 'lat', 'lon');
Ishii_FirstOrder = squeeze(TSLA_AllOrders(:, :, :, 1));
Ishii_lat = lat;
Ishii_lon = lon;
Ishii_mean = squeeze(mean(Ishii_FirstOrder, 3));
Ishii_mean(isnan(Ishii_mean)) = 0;

% 绘制空间分布图
fprintf('绘制空间分布图...\n');

% 创建设置文件
settings_file = 'FirstOrder_Spatial_Settings.mat';
save(settings_file, ...
    'EN4_mean', 'EN4_lat', 'EN4_lon', ...
    'IAP_mean', 'IAP_lat', 'IAP_lon', ...
    'Ishii_mean', 'Ishii_lat', 'Ishii_lon');

fprintf('空间分布数据已保存到: %s\n', settings_file);

% 打印统计信息
fprintf('\n=== 一阶项统计信息 ===\n');
fprintf('EN4: 最小值 = %.4f, 最大值 = %.4f\n', min(EN4_mean(:)), max(EN4_mean(:)));
fprintf('IAP: 最小值 = %.4f, 最大值 = %.4f\n', min(IAP_mean(:)), max(IAP_mean(:)));
fprintf('Ishii: 最小值 = %.4f, 最大值 = %.4f\n', min(Ishii_mean(:)), max(Ishii_mean(:)));

% 找出最大值和最小值的位置
[EN4_max_val, EN4_max_idx] = max(EN4_mean(:));
[EN4_max_lon, EN4_max_lat] = ind2sub(size(EN4_mean), EN4_max_idx);
[EN4_min_val, EN4_min_idx] = min(EN4_mean(:));
[EN4_min_lon, EN4_min_lat] = ind2sub(size(EN4_mean), EN4_min_idx);

[IAP_max_val, IAP_max_idx] = max(IAP_mean(:));
[IAP_max_lon, IAP_max_lat] = ind2sub(size(IAP_mean), IAP_max_idx);
[IAP_min_val, IAP_min_idx] = min(IAP_mean(:));
[IAP_min_lon, IAP_min_lat] = ind2sub(size(IAP_mean), IAP_min_idx);

[Ishii_max_val, Ishii_max_idx] = max(Ishii_mean(:));
[Ishii_max_lon, Ishii_max_lat] = ind2sub(size(Ishii_mean), Ishii_max_idx);
[Ishii_min_val, Ishii_min_idx] = min(Ishii_mean(:));
[Ishii_min_lon, Ishii_min_lat] = ind2sub(size(Ishii_mean), Ishii_min_idx);

fprintf('\n=== 最大值和最小值位置 ===\n');
fprintf('EN4 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', EN4_max_val, EN4_lon(EN4_max_lon), EN4_lat(EN4_max_lat));
fprintf('EN4 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', EN4_min_val, EN4_lon(EN4_min_lon), EN4_lat(EN4_min_lat));
fprintf('IAP 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', IAP_max_val, IAP_lon(IAP_max_lon), IAP_lat(IAP_max_lat));
fprintf('IAP 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', IAP_min_val, IAP_lon(IAP_min_lon), IAP_lat(IAP_min_lat));
fprintf('Ishii 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', Ishii_max_val, Ishii_lon(Ishii_max_lon), Ishii_lat(Ishii_max_lat));
fprintf('Ishii 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', Ishii_min_val, Ishii_lon(Ishii_min_lon), Ishii_lat(Ishii_min_lat));

fprintf('\n>>> 分析完成！<<<\n');
