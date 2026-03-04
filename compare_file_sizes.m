%% 对比新旧脚本的输出文件大小和内容
clear; clc;

% 新脚本输出目录
NewDir = 'D:\work\IAP_mat_data';
% 旧脚本输出目录
OldDir = 'D:\work\MAT_Data';

fprintf('=== 文件大小对比 ===\n\n');

% 对比 T 项
fprintf('【T项 (TSLA)】\n');
NewT = fullfile(NewDir, 'IAP_TSLA_Terms_1to8_Average.mat');
OldT = fullfile(OldDir, 'IAP_TSLA_Terms_1to8_Average.mat');

if exist(NewT, 'file')
    info = dir(NewT);
    fprintf('  新脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end
if exist(OldT, 'file')
    info = dir(OldT);
    fprintf('  旧脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end

% 对比 S 项
fprintf('\n【S项 (SSLA)】\n');
NewS = fullfile(NewDir, 'IAP_S_Terms_1to8_Average.mat');
OldS = fullfile(OldDir, 'IAP_S_Terms_1to8_Average.mat');

if exist(NewS, 'file')
    info = dir(NewS);
    fprintf('  新脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end
if exist(OldS, 'file')
    info = dir(OldS);
    fprintf('  旧脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end

% 对比混合项
fprintf('\n【混合项 (Cross)】\n');
NewC = fullfile(NewDir, 'IAP_Cross_Terms_1to8_Average.mat');
OldC = fullfile(OldDir, 'IAP_Cross_Terms_1to8_Average.mat');

if exist(NewC, 'file')
    info = dir(NewC);
    fprintf('  新脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end
if exist(OldC, 'file')
    info = dir(OldC);
    fprintf('  旧脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end

% 对比标准态版本
fprintf('\n【标准态 T项】\n');
NewStdT = fullfile(NewDir, 'IAP_TSLA_Terms_1to8_StdRef.mat');
OldStdT = fullfile(OldDir, 'IAP_TSLA_Terms_1to8_StdRef.mat');

if exist(NewStdT, 'file')
    info = dir(NewStdT);
    fprintf('  新脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end
if exist(OldStdT, 'file')
    info = dir(OldStdT);
    fprintf('  旧脚本: %.2f MB\n', info.bytes / 1024 / 1024);
end

% 详细分析数据内容
fprintf('\n=== 数据内容分析 ===\n');

if exist(NewT, 'file') && exist(OldT, 'file')
    fprintf('\n【T项数据对比】\n');
    load(NewT, 'TSLA_AllOrders');
    NewData = TSLA_AllOrders;
    load(OldT, 'TSLA_AllOrders', 'time_vec');
    OldData = TSLA_AllOrders;
    
    fprintf('  新脚本维度: %s\n', mat2str(size(NewData)));
    fprintf('  旧脚本维度: %s\n', mat2str(size(OldData)));
    
    fprintf('  新脚本数据类型: %s\n', class(NewData));
    fprintf('  旧脚本数据类型: %s\n', class(OldData));
    
    fprintf('  新脚本NaN比例: %.2f%%\n', 100 * sum(isnan(NewData(:))) / numel(NewData));
    fprintf('  旧脚本NaN比例: %.2f%%\n', 100 * sum(isnan(OldData(:))) / numel(OldData));
    
    % 检查非NaN值的范围
    NewValid = NewData(~isnan(NewData));
    OldValid = OldData(~isnan(OldData));
    
    fprintf('  新脚本非NaN值范围: [%.4f, %.4f]\n', min(NewValid), max(NewValid));
    fprintf('  旧脚本非NaN值范围: [%.4f, %.4f]\n', min(OldValid), max(OldValid));
end

if exist(NewS, 'file') && exist(OldS, 'file')
    fprintf('\n【S项数据对比】\n');
    load(NewS, 'SSLA_AllOrders');
    NewData = SSLA_AllOrders;
    load(OldS, 'SSLA_AllOrders');
    OldData = SSLA_AllOrders;
    
    fprintf('  新脚本维度: %s\n', mat2str(size(NewData)));
    fprintf('  旧脚本维度: %s\n', mat2str(size(OldData)));
    
    fprintf('  新脚本NaN比例: %.2f%%\n', 100 * sum(isnan(NewData(:))) / numel(NewData));
    fprintf('  旧脚本NaN比例: %.2f%%\n', 100 * sum(isnan(OldData(:))) / numel(OldData));
end

fprintf('\n=== 结论 ===\n');
