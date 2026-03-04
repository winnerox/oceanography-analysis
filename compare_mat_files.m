%% 对比两个 mat 文件的结构
clear; clc;

%% 1. 加载两个文件
file1 = 'D:\work\IAP_mat_data\IAP_TSLA_Terms_1to8_Average.mat';
file2 = 'D:\work\IAP_TSLA_Terms\IAP_TSLA_Terms_1to8_Average.mat';

fprintf('========================================\n');
fprintf('对比两个 mat 文件的结构\n');
fprintf('========================================\n\n');

%% 2. 分析文件1
fprintf('文件1: %s\n', file1);
fprintf('文件大小: ');
try
    file_info1 = dir(file1);
    fprintf('%.2f MB\n', file_info1.bytes / 1024 / 1024);
catch
    fprintf('文件不存在\n');
end

fprintf('\n变量列表:\n');
try
    vars1 = whos('-file', file1);
    for i = 1:length(vars1)
        fprintf('  %s: %s, 大小: %s\n', vars1(i).name, vars1(i).class, mat2str(vars1(i).size));
    end
    
    % 加载并显示详细维度
    load(file1, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis');
    fprintf('\n详细维度:\n');
    fprintf('  TSLA_AllOrders: %s\n', mat2str(size(TSLA_AllOrders)));
    fprintf('  lat: %d 点 (%.1f 到 %.1f)\n', length(lat), min(lat), max(lat));
    fprintf('  lon: %d 点 (%.1f 到 %.1f)\n', length(lon), min(lon), max(lon));
    fprintf('  Time_Axis: %d 个时间点\n', length(Time_Axis));
    
    % 保存用于对比
    dim1 = size(TSLA_AllOrders);
    lat1 = lat; lon1 = lon; time1 = Time_Axis;
    
catch ME
    fprintf('加载失败: %s\n', ME.message);
end

%% 3. 分析文件2
fprintf('\n\n========================================\n');
fprintf('文件2: %s\n', file2);
fprintf('文件大小: ');
try
    file_info2 = dir(file2);
    fprintf('%.2f MB\n', file_info2.bytes / 1024 / 1024);
catch
    fprintf('文件不存在\n');
end

fprintf('\n变量列表:\n');
try
    vars2 = whos('-file', file2);
    for i = 1:length(vars2)
        fprintf('  %s: %s, 大小: %s\n', vars2(i).name, vars2(i).class, mat2str(vars2(i).size));
    end
    
    % 加载并显示详细维度
    load(file2, 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis');
    fprintf('\n详细维度:\n');
    fprintf('  TSLA_AllOrders: %s\n', mat2str(size(TSLA_AllOrders)));
    fprintf('  lat: %d 点 (%.1f 到 %.1f)\n', length(lat), min(lat), max(lat));
    fprintf('  lon: %d 点 (%.1f 到 %.1f)\n', length(lon), min(lon), max(lon));
    fprintf('  Time_Axis: %d 个时间点\n', length(Time_Axis));
    
    % 保存用于对比
    dim2 = size(TSLA_AllOrders);
    lat2 = lat; lon2 = lon; time2 = Time_Axis;
    
catch ME
    fprintf('加载失败: %s\n', ME.message);
end

%% 4. 对比差异
fprintf('\n\n========================================\n');
fprintf('差异对比\n');
fprintf('========================================\n');

if exist('dim1', 'var') && exist('dim2', 'var')
    fprintf('\nTSLA_AllOrders 维度对比:\n');
    fprintf('  文件1: %s\n', mat2str(dim1));
    fprintf('  文件2: %s\n', mat2str(dim2));
    fprintf('  差异: %s\n', mat2str(dim1 - dim2));
    
    fprintf('\n坐标对比:\n');
    fprintf('  lat 长度: %d vs %d\n', length(lat1), length(lat2));
    fprintf('  lon 长度: %d vs %d\n', length(lon1), length(lon2));
    fprintf('  Time_Axis 长度: %d vs %d\n', length(time1), length(time2));
    
    % 检查坐标是否相同
    if length(lat1) == length(lat2) && isequal(lat1, lat2)
        fprintf('  ✅ lat 坐标相同\n');
    else
        fprintf('  ⚠️ lat 坐标不同\n');
    end
    
    if length(lon1) == length(lon2) && isequal(lon1, lon2)
        fprintf('  ✅ lon 坐标相同\n');
    else
        fprintf('  ⚠️ lon 坐标不同\n');
    end
    
    if length(time1) == length(time2) && isequal(time1, time2)
        fprintf('  ✅ Time_Axis 相同\n');
    else
        fprintf('  ⚠️ Time_Axis 不同\n');
    end
end

fprintf('\n✅ 对比完成！\n');