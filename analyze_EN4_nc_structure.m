%% EN4 NetCDF 文件结构分析脚本
clear; clc;

%% 1. 查找 EN4 样本文件
fprintf('========================================\n');
fprintf('查找 EN4 样本文件:\n');
fprintf('========================================\n');

DataRoot = 'D:\work\EN4_analyses_c13_last20years';
Years = 2005:2024;

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

if isempty(SampleFile)
    error('❌ 未找到 EN4 样本文件');
end

fprintf('样本文件: %s\n', SampleFile);

%% 2. 分析 EN4 文件结构
fprintf('\n========================================\n');
fprintf('分析 EN4 文件结构:\n');
fprintf('========================================\n');

if ~exist(SampleFile, 'file')
    error('❌ 文件不存在: %s', SampleFile);
end

% 获取文件信息
file_info = ncinfo(SampleFile);
fprintf('文件: %s\n', SampleFile);
fprintf('维度数量: %d\n', length(file_info.Dimensions));
fprintf('变量数量: %d\n', length(file_info.Variables));

% 显示维度信息
fprintf('\n维度信息:\n');
for i = 1:length(file_info.Dimensions)
    dim = file_info.Dimensions(i);
    fprintf('  %s: %d\n', dim.Name, dim.Length);
end

% 显示变量信息
fprintf('\n变量信息:\n');
for i = 1:length(file_info.Variables)
    var = file_info.Variables(i);
    fprintf('  %s: ', var.Name);
    if ~isempty(var.Dimensions)
        fprintf('维度: ');
        for j = 1:length(var.Dimensions)
            fprintf('%s(%d) ', var.Dimensions(j).Name, var.Dimensions(j).Length);
        end
    end
    fprintf('大小: %s\n', mat2str(var.Size));
end

%% 3. 读取并显示实际数据维度
fprintf('\n========================================\n');
fprintf('实际数据读取测试:\n');
fprintf('========================================\n');

try
    % 读取坐标变量
    lon = ncread(SampleFile, 'lon');
    lat = ncread(SampleFile, 'lat');
    depth = ncread(SampleFile, 'depth');
    
    fprintf('lon 维度: %d (%.1f 到 %.1f)\n', length(lon), min(lon), max(lon));
    fprintf('lat 维度: %d (%.1f 到 %.1f)\n', length(lat), min(lat), max(lat));
    fprintf('depth 维度: %d (%.1f 到 %.1f m)\n', length(depth), min(depth), max(depth));
    
    % 读取温度和盐度数据
    temp_data = ncread(SampleFile, 'temperature');
    salt_data = ncread(SampleFile, 'salinity');
    
    fprintf('\ntemperature 变量实际维度: %s\n', mat2str(size(temp_data)));
    fprintf('salinity 变量实际维度: %s\n', mat2str(size(salt_data)));
    
    % 显示数据范围
    fprintf('\ntemperature 数据范围: %.2f 到 %.2f K\n', min(temp_data(:)), max(temp_data(:)));
    fprintf('salinity 数据范围: %.2f 到 %.2f PSU\n', min(salt_data(:)), max(salt_data(:)));
    
    % 检测温度单位 (开尔文 vs 摄氏度)
    temp_mean = mean(temp_data(:), 'omitnan');
    if temp_mean > 100
        fprintf('\n⚠️ 温度单位检测: 开尔文 (K)\n');
        fprintf('   平均值: %.2f K = %.2f °C\n', temp_mean, temp_mean - 273.15);
    else
        fprintf('\n✅ 温度单位检测: 摄氏度 (°C)\n');
        fprintf('   平均值: %.2f °C\n', temp_mean);
    end
    
catch ME
    fprintf('读取错误: %s\n', ME.message);
end

%% 4. 分析维度顺序
fprintf('\n========================================\n');
fprintf('维度顺序分析:\n');
fprintf('========================================\n');

if exist('temp_data', 'var') && exist('lon', 'var') && exist('lat', 'var') && exist('depth', 'var')
    data_sizes = size(temp_data);
    fprintf('原始数据维度: %s\n', mat2str(data_sizes));
    fprintf('lon 长度: %d\n', length(lon));
    fprintf('lat 长度: %d\n', length(lat));
    fprintf('depth 长度: %d\n', length(depth));
    
    % 推断维度顺序
    fprintf('\n推断维度顺序:\n');
    
    % 匹配各维度
    dim_match = zeros(1, 3);
    for i = 1:3
        if data_sizes(i) == length(lon)
            dim_match(i) = 1;  % lon
        elseif data_sizes(i) == length(lat)
            dim_match(i) = 2;  % lat
        elseif data_sizes(i) == length(depth)
            dim_match(i) = 3;  % depth
        end
    end
    
    dim_names = {'lon','lat','depth'};
    fprintf('  维度匹配: [%s]\n', strjoin(dim_names(dim_match), ', '));
    
    if all(dim_match == [1, 2, 3])
        fprintf('  推断: [lon, lat, depth]\n');
    elseif all(dim_match == [2, 1, 3])
        fprintf('  推断: [lat, lon, depth]\n');
    elseif all(dim_match == [1, 3, 2])
        fprintf('  推断: [lon, depth, lat]\n');
    elseif all(dim_match == [3, 1, 2])
        fprintf('  推断: [depth, lon, lat]\n');
    elseif all(dim_match == [2, 3, 1])
        fprintf('  推断: [lat, depth, lon]\n');
    elseif all(dim_match == [3, 2, 1])
        fprintf('  推断: [depth, lat, lon]\n');
    else
        fprintf('  ⚠️ 无法确定维度顺序\n');
    end
end

%% 5. 测试不同的 permute 组合
fprintf('\n========================================\n');
fprintf('测试 permute 组合:\n');
fprintf('========================================\n');

if exist('temp_data', 'var')
    test_sizes = size(temp_data);
    fprintf('原始维度: %s\n', mat2str(test_sizes));
    
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
end

%% 6. 建议的读取方式
fprintf('\n========================================\n');
fprintf('建议的读取方式:\n');
fprintf('========================================\n');

if exist('lon', 'var') && exist('lat', 'var') && exist('depth', 'var')
    fprintf('根据分析，建议使用以下读取方式:\n');
    fprintf('1. 原始数据维度: [lon, lat, depth, time]\n');
    fprintf('2. 目标处理维度: [lat, depth, lon] (便于深度积分)\n');
    fprintf('3. 因此 permute 顺序应为: [2, 3, 1]\n');
    fprintf('\n代码示例:\n');
    fprintf('  Nz = length(depth);\n');
    fprintf('  temp_data = ncread(file, ''temperature'', [1 1 1 1], [Inf Inf Nz 1]);\n');
    fprintf('  %% 原始维度 [lon, lat, depth] → 目标 [lat, depth, lon]\n');
    fprintf('  temp_permuted = permute(temp_data, [2, 3, 1]);\n');
    fprintf('  size(temp_permuted) = %s\n', mat2str([length(lat), length(depth), length(lon)]));
    
    fprintf('\n⚠️ EN4 特有注意事项:\n');
    fprintf('  - 温度单位为开尔文 (K)，需要转换为摄氏度: T_C = T_K - 273.15\n');
    fprintf('  - 文件包含时间维度 (第4维)，读取时需指定 [1 1 1 1] 起始位置\n');
end

fprintf('\n✅ EN4 文件结构分析完成！\n');
