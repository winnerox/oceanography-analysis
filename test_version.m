%% 测试当前 test.m 版本
clear; clc;

% 读取 test.m 文件
fid = fopen('test.m', 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = lines{1};

% 查找关键行
for i = 1:min(100, length(lines))
    if contains(lines{i}, 'for k = 2:Nz-1')
        fprintf('找到旧版本代码在第 %d 行: %s\n', i, lines{i});
    elseif contains(lines{i}, 'for k = 2:num_depth-1')
        fprintf('找到新版本代码在第 %d 行: %s\n', i, lines{i});
    end
end

% 检查第59行内容
if length(lines) >= 59
    fprintf('\n第59行内容: %s\n', lines{59});
else
    fprintf('\n文件少于59行\n');
end

% 检查第62-66行内容
fprintf('\n检查积分权重计算部分:\n');
for i = 60:70
    if i <= length(lines)
        fprintf('%3d: %s\n', i, lines{i});
    end
end