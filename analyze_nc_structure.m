%% NetCDF 文件结构分析脚本
clear; clc;

%% 1. 分析温度文件
fprintf('========================================\n');
fprintf('分析温度文件:\n');
fprintf('========================================\n');
temp_file = 'D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc';

if ~exist(temp_file, 'file')
    error('❌ 温度文件不存在: %s', temp_file);
end

% 获取文件信息
temp_info = ncinfo(temp_file);
fprintf('文件: %s\n', temp_file);
fprintf('维度数量: %d\n', length(temp_info.Dimensions));
fprintf('变量数量: %d\n', length(temp_info.Variables));

% 显示维度信息
fprintf('\n维度信息:\n');
for i = 1:length(temp_info.Dimensions)
    dim = temp_info.Dimensions(i);
    fprintf('  %s: %d\n', dim.Name, dim.Length);
end

% 显示变量信息
fprintf('\n变量信息:\n');
for i = 1:length(temp_info.Variables)
    var = temp_info.Variables(i);
    fprintf('  %s: ', var.Name);
    if ~isempty(var.Dimensions)
        fprintf('维度: ');
        for j = 1:length(var.Dimensions)
            fprintf('%s(%d) ', var.Dimensions(j).Name, var.Dimensions(j).Length);
        end
    end
    fprintf('大小: %s\n', mat2str(var.Size));
end

% 读取并显示实际数据维度
fprintf('\n实际数据读取测试:\n');
try
    temp_data = ncread(temp_file, 'temp');
    fprintf('  temp 变量实际维度: %s\n', mat2str(size(temp_data)));
    
    % 读取坐标变量
    lon_temp = ncread(temp_file, 'lon');
    lat_temp = ncread(temp_file, 'lat');
    depth_temp = ncread(temp_file, 'depth_std');
    
    fprintf('  lon 维度: %d (%.1f 到 %.1f)\n', length(lon_temp), min(lon_temp), max(lon_temp));
    fprintf('  lat 维度: %d (%.1f 到 %.1f)\n', length(lat_temp), min(lat_temp), max(lat_temp));
    fprintf('  depth_std 维度: %d (%.1f 到 %.1f)\n', length(depth_temp), min(depth_temp), max(depth_temp));
    
    % 显示数据范围
    fprintf('  temp 数据范围: %.2f 到 %.2f\n', min(temp_data(:)), max(temp_data(:)));
    
catch ME
    fprintf('  读取错误: %s\n', ME.message);
end

%% 2. 分析盐度文件
fprintf('\n\n========================================\n');
fprintf('分析盐度文件:\n');
fprintf('========================================\n');
salt_file = 'D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc';

if ~exist(salt_file, 'file')
    error('❌ 盐度文件不存在: %s', salt_file);
end

% 获取文件信息
salt_info = ncinfo(salt_file);
fprintf('文件: %s\n', salt_file);
fprintf('维度数量: %d\n', length(salt_info.Dimensions));
fprintf('变量数量: %d\n', length(salt_info.Variables));

% 显示维度信息
fprintf('\n维度信息:\n');
for i = 1:length(salt_info.Dimensions)
    dim = salt_info.Dimensions(i);
    fprintf('  %s: %d\n', dim.Name, dim.Length);
end

% 显示变量信息
fprintf('\n变量信息:\n');
for i = 1:length(salt_info.Variables)
    var = salt_info.Variables(i);
    fprintf('  %s: ', var.Name);
    if ~isempty(var.Dimensions)
        fprintf('维度: ');
        for j = 1:length(var.Dimensions)
            fprintf('%s(%d) ', var.Dimensions(j).Name, var.Dimensions(j).Length);
        end
    end
    fprintf('大小: %s\n', mat2str(var.Size));
end

% 读取并显示实际数据维度
fprintf('\n实际数据读取测试:\n');
try
    salt_data = ncread(salt_file, 'salinity');
    fprintf('  salinity 变量实际维度: %s\n', mat2str(size(salt_data)));
    
    % 读取坐标变量
    lon_salt = ncread(salt_file, 'lon');
    lat_salt = ncread(salt_file, 'lat');
    depth_salt = ncread(salt_file, 'depth_std');
    
    fprintf('  lon 维度: %d (%.1f 到 %.1f)\n', length(lon_salt), min(lon_salt), max(lon_salt));
    fprintf('  lat 维度: %d (%.1f 到 %.1f)\n', length(lat_salt), min(lat_salt), max(lat_salt));
    fprintf('  depth_std 维度: %d (%.1f 到 %.1f)\n', length(depth_salt), min(depth_salt), max(depth_salt));
    
    % 显示数据范围
    fprintf('  salinity 数据范围: %.2f 到 %.2f\n', min(salt_data(:)), max(salt_data(:)));
    
catch ME
    fprintf('  读取错误: %s\n', ME.message);
end

%% 3. 维度一致性检查
fprintf('\n\n========================================\n');
fprintf('维度一致性检查:\n');
fprintf('========================================\n');

% 检查坐标是否一致
lon_match = isequal(lon_temp, lon_salt);
lat_match = isequal(lat_temp, lat_salt);
depth_match = isequal(depth_temp, depth_salt);

fprintf('lon 坐标一致: %s\n', mat2str(lon_match));
fprintf('lat 坐标一致: %s\n', mat2str(lat_match));
fprintf('depth 坐标一致: %s\n', mat2str(depth_match));

if lon_match && lat_match && depth_match
    fprintf('✅ 所有坐标维度完全一致！\n');
else
    fprintf('⚠️ 坐标维度不一致，需要调整！\n');
end

% 检查数据维度顺序
fprintf('\n数据维度顺序分析:\n');
fprintf('温度数据维度: %s\n', mat2str(size(temp_data)));
fprintf('盐度数据维度: %s\n', mat2str(size(salt_data)));

% 分析维度顺序
if length(size(temp_data)) == 3
    fprintf('\n推断维度顺序:\n');
    fprintf('  根据变量名和常见约定:\n');
    fprintf('  - lon: 经度 (通常第一维)\n');
    fprintf('  - lat: 纬度 (通常第二维)\n');
    fprintf('  - depth_std: 深度 (通常第三维)\n');
    fprintf('  所以数据维度顺序可能是: [lon, lat, depth]\n');
end

%% 4. 测试不同的 permute 组合
fprintf('\n\n========================================\n');
fprintf('测试 permute 组合:\n');
fprintf('========================================\n');

test_sizes = size(temp_data);
fprintf('原始维度: %s\n', mat2str(test_sizes));

% 测试不同的 permute 组合
permute_combinations = {
    [1, 2, 3], '原始顺序';
    [2, 1, 3], '交换 lon 和 lat';
    [3, 1, 2], '深度第一维';
    [3, 2, 1], '深度第一维，交换 lon 和 lat';
    [1, 3, 2], '交换 lat 和 depth';
    [2, 3, 1], '交换所有维度';
};

for i = 1:size(permute_combinations, 1)
    perm_order = permute_combinations{i, 1};
    description = permute_combinations{i, 2};
    
    test_permuted = permute(temp_data(1:min(2, test_sizes(1)), ...
                                      1:min(2, test_sizes(2)), ...
                                      1:min(2, test_sizes(3))), perm_order);
    
    fprintf('  permute(%s): %s → 维度: %s\n', ...
        mat2str(perm_order), description, mat2str(size(test_permuted)));
end

%% 5. 建议的读取方式
fprintf('\n\n========================================\n');
fprintf('建议的读取方式:\n');
fprintf('========================================\n');

fprintf('根据分析，建议使用以下读取方式:\n');
fprintf('1. 原始数据维度: [lon, lat, depth]\n');
fprintf('2. 目标处理维度: [lat, depth, lon] (便于深度积分)\n');
fprintf('3. 因此 permute 顺序应为: [2, 3, 1]\n');
fprintf('\n代码示例:\n');
fprintf('  temp_data = ncread(temp_file, ''temp'', [1 1 1], [Inf Inf Inf]);\n');
fprintf('  temp_permuted = permute(temp_data, [2, 3, 1]); %% [lat, depth, lon]\n');
fprintf('  size(temp_permuted) = %s\n', mat2str([length(lat_temp), length(depth_temp), length(lon_temp)]));

fprintf('\n✅ 分析完成！\n');