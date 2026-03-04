%% Analyze_FirstOrder_Spatial.m
% =========================================================================
% 功能：分析三个数据集（EN4、IAP、Ishii）的平均态一阶项空间分布
% 输出：详细的空间分布分析报告
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

% 计算统计信息
fprintf('\n=== 一阶项统计信息 ===\n');
fprintf('EN4: 最小值 = %.4f, 最大值 = %.4f\n', min(EN4_mean(:)), max(EN4_mean(:)));
fprintf('IAP: 最小值 = %.4f, 最大值 = %.4f\n', min(IAP_mean(:)), max(IAP_mean(:)));
fprintf('Ishii: 最小值 = %.4f, 最大值 = %.4f\n', min(Ishii_mean(:)), max(Ishii_mean(:)));

% 找出最大值和最小值的位置
fprintf('\n=== 最大值和最小值位置 ===\n');

% EN4
[EN4_max_val, EN4_max_idx] = max(EN4_mean(:));
[EN4_max_lat_idx, EN4_max_lon_idx] = ind2sub(size(EN4_mean), EN4_max_idx);
% 确保索引在有效范围内
EN4_max_lon_idx = min(EN4_max_lon_idx, length(EN4_lon));
EN4_max_lat_idx = min(EN4_max_lat_idx, length(EN4_lat));
EN4_max_lon = EN4_lon(EN4_max_lon_idx);
EN4_max_lat = EN4_lat(EN4_max_lat_idx);
fprintf('EN4 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', EN4_max_val, EN4_max_lon, EN4_max_lat);

[EN4_min_val, EN4_min_idx] = min(EN4_mean(:));
[EN4_min_lat_idx, EN4_min_lon_idx] = ind2sub(size(EN4_mean), EN4_min_idx);
% 确保索引在有效范围内
EN4_min_lon_idx = min(EN4_min_lon_idx, length(EN4_lon));
EN4_min_lat_idx = min(EN4_min_lat_idx, length(EN4_lat));
EN4_min_lon = EN4_lon(EN4_min_lon_idx);
EN4_min_lat = EN4_lat(EN4_min_lat_idx);
fprintf('EN4 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', EN4_min_val, EN4_min_lon, EN4_min_lat);

% IAP
[IAP_max_val, IAP_max_idx] = max(IAP_mean(:));
[IAP_max_lat_idx, IAP_max_lon_idx] = ind2sub(size(IAP_mean), IAP_max_idx);
% 确保索引在有效范围内
IAP_max_lon_idx = min(IAP_max_lon_idx, length(IAP_lon));
IAP_max_lat_idx = min(IAP_max_lat_idx, length(IAP_lat));
IAP_max_lon = IAP_lon(IAP_max_lon_idx);
IAP_max_lat = IAP_lat(IAP_max_lat_idx);
fprintf('IAP 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', IAP_max_val, IAP_max_lon, IAP_max_lat);

[IAP_min_val, IAP_min_idx] = min(IAP_mean(:));
[IAP_min_lat_idx, IAP_min_lon_idx] = ind2sub(size(IAP_mean), IAP_min_idx);
% 确保索引在有效范围内
IAP_min_lon_idx = min(IAP_min_lon_idx, length(IAP_lon));
IAP_min_lat_idx = min(IAP_min_lat_idx, length(IAP_lat));
IAP_min_lon = IAP_lon(IAP_min_lon_idx);
IAP_min_lat = IAP_lat(IAP_min_lat_idx);
fprintf('IAP 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', IAP_min_val, IAP_min_lon, IAP_min_lat);

% Ishii
[Ishii_max_val, Ishii_max_idx] = max(Ishii_mean(:));
[Ishii_max_lat_idx, Ishii_max_lon_idx] = ind2sub(size(Ishii_mean), Ishii_max_idx);
% 确保索引在有效范围内
Ishii_max_lon_idx = min(Ishii_max_lon_idx, length(Ishii_lon));
Ishii_max_lat_idx = min(Ishii_max_lat_idx, length(Ishii_lat));
Ishii_max_lon = Ishii_lon(Ishii_max_lon_idx);
Ishii_max_lat = Ishii_lat(Ishii_max_lat_idx);
fprintf('Ishii 最大值: %.4f mm at lon=%.2f, lat=%.2f\n', Ishii_max_val, Ishii_max_lon, Ishii_max_lat);

[Ishii_min_val, Ishii_min_idx] = min(Ishii_mean(:));
[Ishii_min_lat_idx, Ishii_min_lon_idx] = ind2sub(size(Ishii_mean), Ishii_min_idx);
% 确保索引在有效范围内
Ishii_min_lon_idx = min(Ishii_min_lon_idx, length(Ishii_lon));
Ishii_min_lat_idx = min(Ishii_min_lat_idx, length(Ishii_lat));
Ishii_min_lon = Ishii_lon(Ishii_min_lon_idx);
Ishii_min_lat = Ishii_lat(Ishii_min_lat_idx);
fprintf('Ishii 最小值: %.4f mm at lon=%.2f, lat=%.2f\n', Ishii_min_val, Ishii_min_lon, Ishii_min_lat);

% 找出值特别大的区域（超过最大值的90%）
fprintf('\n=== 值特别大的区域 ===\n');

% EN4
EN4_threshold = max(EN4_mean(:)) * 0.9;
EN4_large_areas = find(EN4_mean > EN4_threshold);
fprintf('EN4: 有 %d 个点超过最大值的90%% (%.4f mm)\n', length(EN4_large_areas), EN4_threshold);
if length(EN4_large_areas) > 0
    for i = 1:min(5, length(EN4_large_areas))
        [lat_idx, lon_idx] = ind2sub(size(EN4_mean), EN4_large_areas(i));
        fprintf('  - lon=%.2f, lat=%.2f: %.4f mm\n', EN4_lon(lon_idx), EN4_lat(lat_idx), EN4_mean(lat_idx, lon_idx));
    end
end

% IAP
IAP_threshold = max(IAP_mean(:)) * 0.9;
IAP_large_areas = find(IAP_mean > IAP_threshold);
fprintf('IAP: 有 %d 个点超过最大值的90%% (%.4f mm)\n', length(IAP_large_areas), IAP_threshold);
if length(IAP_large_areas) > 0
    for i = 1:min(5, length(IAP_large_areas))
        [lat_idx, lon_idx] = ind2sub(size(IAP_mean), IAP_large_areas(i));
        fprintf('  - lon=%.2f, lat=%.2f: %.4f mm\n', IAP_lon(lon_idx), IAP_lat(lat_idx), IAP_mean(lat_idx, lon_idx));
    end
end

% Ishii
Ishii_threshold = max(Ishii_mean(:)) * 0.9;
Ishii_large_areas = find(Ishii_mean > Ishii_threshold);
fprintf('Ishii: 有 %d 个点超过最大值的90%% (%.4f mm)\n', length(Ishii_large_areas), Ishii_threshold);
if length(Ishii_large_areas) > 0
    for i = 1:min(5, length(Ishii_large_areas))
        [lat_idx, lon_idx] = ind2sub(size(Ishii_mean), Ishii_large_areas(i));
        fprintf('  - lon=%.2f, lat=%.2f: %.4f mm\n', Ishii_lon(lon_idx), Ishii_lat(lat_idx), Ishii_mean(lat_idx, lon_idx));
    end
end

% 绘制3D空间分布图
fprintf('\n绘制3D空间分布图...\n');

% 检查并调整数据维度
fprintf('检查数据维度...\n');
fprintf('EN4_mean size: %s\n', mat2str(size(EN4_mean)));
fprintf('EN4_lon length: %d\n', length(EN4_lon));
fprintf('EN4_lat length: %d\n', length(EN4_lat));
fprintf('IAP_mean size: %s\n', mat2str(size(IAP_mean)));
fprintf('IAP_lon length: %d\n', length(IAP_lon));
fprintf('IAP_lat length: %d\n', length(IAP_lat));
fprintf('Ishii_mean size: %s\n', mat2str(size(Ishii_mean)));
fprintf('Ishii_lon length: %d\n', length(Ishii_lon));
fprintf('Ishii_lat length: %d\n', length(Ishii_lat));

% 确保维度匹配并绘制3D图
% EN4 3D图
figure('Name', 'EN4 一阶项空间分布 (3D)', 'Position', [100, 100, 1000, 800]);
% 转置数据以匹配meshgrid维度
EN4_mean_transposed = EN4_mean';
[X, Y] = meshgrid(EN4_lon, EN4_lat);
surf(X, Y, EN4_mean_transposed);
title('EN4 一阶项空间分布 (3D)');
xlabel('经度');
ylabel('纬度');
zlabel('值 (mm)');
colorbar;
view(30, 45);
fprintf('EN4 3D图绘制完成\n');

% IAP 3D图
figure('Name', 'IAP 一阶项空间分布 (3D)', 'Position', [100, 100, 1000, 800]);
% 转置数据以匹配meshgrid维度
IAP_mean_transposed = IAP_mean';
[X, Y] = meshgrid(IAP_lon, IAP_lat);
surf(X, Y, IAP_mean_transposed);
title('IAP 一阶项空间分布 (3D)');
xlabel('经度');
ylabel('纬度');
zlabel('值 (mm)');
colorbar;
view(30, 45);
fprintf('IAP 3D图绘制完成\n');

% Ishii 3D图
figure('Name', 'Ishii 一阶项空间分布 (3D)', 'Position', [100, 100, 1000, 800]);
% 转置数据以匹配meshgrid维度
Ishii_mean_transposed = Ishii_mean';
[X, Y] = meshgrid(Ishii_lon, Ishii_lat);
surf(X, Y, Ishii_mean_transposed);
title('Ishii 一阶项空间分布 (3D)');
xlabel('经度');
ylabel('纬度');
zlabel('值 (mm)');
colorbar;
view(30, 45);
fprintf('Ishii 3D图绘制完成\n');

fprintf('\n>>> 分析完成！<<<\n');