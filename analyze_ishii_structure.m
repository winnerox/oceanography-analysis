%% Ishii 数据文件结构分析
clear; clc;

%% 1. 分析温度文件
fprintf('========================================\n');
fprintf('分析 Ishii 温度文件:\n');
fprintf('========================================\n');
temp_file = 'D:\work\Ishii_05_24\Temperature\temp.2005.nc';

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
    % 尝试读取温度数据
    if any(strcmp({temp_info.Variables.Name}, 'temp'))
        temp_data = ncread(temp_file, 'temp');
        fprintf('  temp 变量实际维度: %s\n', mat2str(size(temp_data)));
        fprintf('  temp 数据范围: %.2f 到 %.2f\n', min(temp_data(:)), max(temp_data(:)));
    elseif any(strcmp({temp_info.Variables.Name}, 'temperature'))
        temp_data = ncread(temp_file, 'temperature');
        fprintf('  temperature 变量实际维度: %s\n', mat2str(size(temp_data)));
        fprintf('  temperature 数据范围: %.2f 到 %.2f\n', min(temp_data(:)), max(temp_data(:)));
    else
        fprintf('  未找到温度变量，可用变量: ');
        for i = 1:length(temp_info.Variables)
            fprintf('%s ', temp_info.Variables(i).Name);
        end
        fprintf('\n');
    end
    
    % 读取坐标变量
    coord_vars = {'lon', 'lat', 'depth', 'lev', 'z', 'time', 'month', 'longitude', 'latitude', 'level'};
    for var_name = coord_vars
        if any(strcmp({temp_info.Variables.Name}, var_name{1}))
            var_data = ncread(temp_file, var_name{1});
            fprintf('  %s 维度: %d (%.1f 到 %.1f)\n', var_name{1}, length(var_data), min(var_data(:)), max(var_data(:)));
        end
    end
    
catch ME
    fprintf('  读取错误: %s\n', ME.message);
end

%% 2. 分析盐度文件
fprintf('\n\n========================================\n');
fprintf('分析 Ishii 盐度文件:\n');
fprintf('========================================\n');
salt_file = 'D:\work\Ishii_05_24\Salinity\sal.2005.nc';

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
    % 尝试读取盐度数据
    if any(strcmp({salt_info.Variables.Name}, 'sal'))
        salt_data = ncread(salt_file, 'sal');
        fprintf('  sal 变量实际维度: %s\n', mat2str(size(salt_data)));
        fprintf('  sal 数据范围: %.2f 到 %.2f\n', min(salt_data(:)), max(salt_data(:)));
    elseif any(strcmp({salt_info.Variables.Name}, 'salinity'))
        salt_data = ncread(salt_file, 'salinity');
        fprintf('  salinity 变量实际维度: %s\n', mat2str(size(salt_data)));
        fprintf('  salinity 数据范围: %.2f 到 %.2f\n', min(salt_data(:)), max(salt_data(:)));
    else
        fprintf('  未找到盐度变量，可用变量: ');
        for i = 1:length(salt_info.Variables)
            fprintf('%s ', salt_info.Variables(i).Name);
        end
        fprintf('\n');
    end
    
    % 读取坐标变量
    coord_vars = {'lon', 'lat', 'depth', 'lev', 'z', 'time', 'month', 'longitude', 'latitude', 'level'};
    for var_name = coord_vars
        if any(strcmp({salt_info.Variables.Name}, var_name{1}))
            var_data = ncread(salt_file, var_name{1});
            fprintf('  %s 维度: %d (%.1f 到 %.1f)\n', var_name{1}, length(var_data), min(var_data(:)), max(var_data(:)));
        end
    end
    
catch ME
    fprintf('  读取错误: %s\n', ME.message);
end

%% 3. 维度一致性检查
fprintf('\n\n========================================\n');
fprintf('维度一致性检查:\n');
fprintf('========================================\n');

try
    % 温度文件的坐标
    temp_lon = ncread(temp_file, 'longitude');
    temp_lat = ncread(temp_file, 'latitude');
    temp_level = ncread(temp_file, 'level');
    
    % 盐度文件的坐标
    salt_lon = ncread(salt_file, 'longitude');
    salt_lat = ncread(salt_file, 'latitude');
    salt_level = ncread(salt_file, 'level');
    
    lon_match = isequal(temp_lon, salt_lon);
    lat_match = isequal(temp_lat, salt_lat);
    level_match = isequal(temp_level, salt_level);
    
    fprintf('longitude 坐标一致: %s\n', mat2str(lon_match));
    fprintf('latitude 坐标一致: %s\n', mat2str(lat_match));
    fprintf('level 坐标一致: %s\n', mat2str(level_match));
    
    if lon_match && lat_match && level_match
        fprintf('✅ 所有坐标维度完全一致！\n');
    else
        fprintf('⚠️ 坐标维度不一致，需要调整！\n');
        
        % 显示差异
        if ~lon_match
            fprintf('  longitude 差异: 温度 %d 点, 盐度 %d 点\n', length(temp_lon), length(salt_lon));
        end
        if ~lat_match
            fprintf('  latitude 差异: 温度 %d 点, 盐度 %d 点\n', length(temp_lat), length(salt_lat));
        end
        if ~level_match
            fprintf('  level 差异: 温度 %d 层, 盐度 %d 层\n', length(temp_level), length(salt_level));
        end
    end
    
    % 显示 level 的具体值
    fprintf('\nlevel 层次值:\n');
    fprintf('  温度文件: %s\n', mat2str(temp_level(1:min(10, length(temp_level)))));
    fprintf('  盐度文件: %s\n', mat2str(salt_level(1:min(10, length(salt_level)))));
    
catch ME
    fprintf('❌ 维度一致性检查失败: %s\n', ME.message);
end

%% 4. 数据维度顺序分析
fprintf('\n\n========================================\n');
fprintf('数据维度顺序分析:\n');
fprintf('========================================\n');

fprintf('Ishii 数据维度顺序: [longitude, latitude, level, time]\n');
fprintf('即: [lon, lat, depth, time] = [360, 180, 28, 12]\n');
fprintf('\n建议的 permute 顺序:\n');
fprintf('  目标维度: [lat, depth, lon, time]\n');
fprintf('  permute(..., [2, 3, 1, 4])\n');
fprintf('  结果维度: [180, 28, 360, 12]\n');

fprintf('\n✅ 分析完成！\n');