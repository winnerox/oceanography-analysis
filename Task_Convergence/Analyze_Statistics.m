%% Analyze_Statistics.m
% =========================================================================
% 功能：分析三套数据集（EN4、IAP、Ishii）的泰勒展开收敛性统计信息
% 输出：各数据集的详细统计分析，包括分位数、极值、均值等
% =========================================================================
clear; clc;
addpath('D:\work');

% 定义数据目录
EN4_DataDir = 'D:\work\EN4_TSLA_Terms';
IAP_DataDir = 'D:\work\IAP_mat_data';
Ishii_DataDir = 'D:\work\Ishii_mat_data';

% 加载收敛性测试结果
fprintf('加载收敛性测试结果...\n');
load(fullfile(EN4_DataDir, 'Taylor_Convergence_Test.mat'), 'EN4_RMS_Avg', 'EN4_RMS_Std', 'EN4_Residual_Avg', 'EN4_Residual_Std');
load(fullfile(IAP_DataDir, 'Taylor_Convergence_Test.mat'), 'IAP_RMS_Avg', 'IAP_RMS_Std', 'IAP_Residual_Avg', 'IAP_Residual_Std');
load(fullfile(Ishii_DataDir, 'Taylor_Convergence_Test.mat'), 'Ishii_RMS_Avg', 'Ishii_RMS_Std', 'Ishii_Residual_Avg', 'Ishii_Residual_Std');

% 加载原始数据（用于更详细的分析）
fprintf('加载原始数据...\n');

% 加载 EN4 平均态数据
load(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
EN4_Exact_Avg = TSLA_Exact_Avg;
load(fullfile(EN4_DataDir, 'EN4_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
EN4_TSLA_Avg = TSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
EN4_SSLA_Avg = SSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
EN4_HSLA_Avg = EN4_TSLA_Avg + EN4_SSLA_Avg + Cross_AllOrders;

% 加载 EN4 标准态数据
load(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Std.mat'), 'TSLA_Exact_Std');
EN4_Exact_Std = TSLA_Exact_Std;
load(fullfile(EN4_DataDir, 'EN4_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders');
EN4_TSLA_Std = TSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders');
EN4_SSLA_Std = SSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders');
EN4_HSLA_Std = EN4_TSLA_Std + EN4_SSLA_Std + Cross_AllOrders;

% 加载 IAP 平均态数据
load(fullfile(IAP_DataDir, 'IAP_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
IAP_Exact_Avg = TSLA_Exact_Avg;
load(fullfile(IAP_DataDir, 'IAP_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
IAP_TSLA_Avg = TSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
IAP_SSLA_Avg = SSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
IAP_HSLA_Avg = IAP_TSLA_Avg + IAP_SSLA_Avg + Cross_AllOrders;

% 加载 IAP 标准态数据
load(fullfile(IAP_DataDir, 'IAP_Formula11_Exact_StdRef.mat'), 'TSLA_Exact_StdRef');
IAP_Exact_Std = TSLA_Exact_StdRef;
load(fullfile(IAP_DataDir, 'IAP_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders');
IAP_TSLA_Std = TSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders');
IAP_SSLA_Std = SSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders');
IAP_HSLA_Std = IAP_TSLA_Std + IAP_SSLA_Std + Cross_AllOrders;

% 加载 Ishii 平均态数据
load(fullfile(Ishii_DataDir, 'Ishii_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
Ishii_Exact_Avg = TSLA_Exact_Avg;
load(fullfile(Ishii_DataDir, 'Ishii_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
Ishii_TSLA_Avg = TSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
Ishii_SSLA_Avg = SSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
Ishii_HSLA_Avg = Ishii_TSLA_Avg + Ishii_SSLA_Avg + Cross_AllOrders;

% 加载 Ishii 标准态数据
load(fullfile(Ishii_DataDir, 'Ishii_Formula11_Exact_StdRef.mat'), 'TSLA_Exact_StdRef');
Ishii_Exact_Std = TSLA_Exact_StdRef;
load(fullfile(Ishii_DataDir, 'Ishii_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders');
Ishii_TSLA_Std = TSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders');
Ishii_SSLA_Std = SSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders');
Ishii_HSLA_Std = Ishii_TSLA_Std + Ishii_SSLA_Std + Cross_AllOrders;

% 分析各数据集的平均态
fprintf('\n========== 三套数据平均态统计分析 ==========\n');
fprintf('\n【子图1：第1-3阶 - SSLA, TSLA, HSLA】\n');
fprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n', ...
    '类别', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');

% EN4 低阶（1-3阶）
fprintf('\n--- EN4 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 低阶（1-3阶）
fprintf('\n--- IAP 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 低阶（1-3阶）
fprintf('\n--- Ishii 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 分析中阶（4-6阶）
fprintf('\n\n【子图2：第4-6阶 - SSLA, TSLA, HSLA】\n');

% EN4 中阶（4-6阶）
fprintf('\n--- EN4 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 中阶（4-6阶）
fprintf('\n--- IAP 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 中阶（4-6阶）
fprintf('\n--- Ishii 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 分析高阶（7-8阶）
fprintf('\n\n【子图3：第7-8阶 - SSLA, TSLA, HSLA】\n');

% EN4 高阶（7-8阶）
fprintf('\n--- EN4 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 高阶（7-8阶）
fprintf('\n--- IAP 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 高阶（7-8阶）
fprintf('\n--- Ishii 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Avg(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 分析各数据集的标准态
fprintf('\n\n========== 三套数据标准态统计分析 ==========\n');
fprintf('\n【子图1：第1-3阶 - SSLA, TSLA, HSLA】\n');
fprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n', ...
    '类别', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');

% EN4 低阶（1-3阶）
fprintf('\n--- EN4 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 低阶（1-3阶）
fprintf('\n--- IAP 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 低阶（1-3阶）
fprintf('\n--- Ishii 低阶 (1-3阶) ---.\n');
for order = 1:3
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 分析中阶（4-6阶）
fprintf('\n\n【子图2：第4-6阶 - SSLA, TSLA, HSLA】\n');

% EN4 中阶（4-6阶）
fprintf('\n--- EN4 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 中阶（4-6阶）
fprintf('\n--- IAP 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 中阶（4-6阶）
fprintf('\n--- Ishii 中阶 (4-6阶) ---.\n');
for order = 4:6
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 分析高阶（7-8阶）
fprintf('\n\n【子图3：第7-8阶 - SSLA, TSLA, HSLA】\n');

% EN4 高阶（7-8阶）
fprintf('\n--- EN4 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('EN4 SSLA %d阶\t', order);
    stats = calculate_statistics(EN4_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 TSLA %d阶\t', order);
    stats = calculate_statistics(EN4_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('EN4 HSLA %d阶\t', order);
    stats = calculate_statistics(EN4_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% IAP 高阶（7-8阶）
fprintf('\n--- IAP 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('IAP SSLA %d阶\t', order);
    stats = calculate_statistics(IAP_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP TSLA %d阶\t', order);
    stats = calculate_statistics(IAP_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('IAP HSLA %d阶\t', order);
    stats = calculate_statistics(IAP_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% Ishii 高阶（7-8阶）
fprintf('\n--- Ishii 高阶 (7-8阶) ---.\n');
for order = 7:8
    fprintf('Ishii SSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_SSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii TSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_TSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
    
    fprintf('Ishii HSLA %d阶\t', order);
    stats = calculate_statistics(Ishii_HSLA_Std(:, :, :, order));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

% 收敛性统计汇总
fprintf('\n\n========== 收敛性统计汇总 ==========\n');
fprintf('\n--- EN4 收敛性 ---.\n');
fprintf('平均态最终残差: %.4f mm (8阶)\n', EN4_RMS_Avg(end));
fprintf('标准态最终残差: %.4f mm (8阶)\n', EN4_RMS_Std(end));
fprintf('平均态收敛比率 (8阶/1阶): %.2f%%\n', (EN4_RMS_Avg(end)/EN4_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (8阶/1阶): %.2f%%\n', (EN4_RMS_Std(end)/EN4_RMS_Std(1))*100);

fprintf('\n--- IAP 收敛性 ---.\n');
fprintf('平均态最终残差: %.4f mm (8阶)\n', IAP_RMS_Avg(end));
fprintf('标准态最终残差: %.4f mm (8阶)\n', IAP_RMS_Std(end));
fprintf('平均态收敛比率 (8阶/1阶): %.2f%%\n', (IAP_RMS_Avg(end)/IAP_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (8阶/1阶): %.2f%%\n', (IAP_RMS_Std(end)/IAP_RMS_Std(1))*100);

fprintf('\n--- Ishii 收敛性 ---.\n');
fprintf('平均态最终残差: %.4f mm (8阶)\n', Ishii_RMS_Avg(end));
fprintf('标准态最终残差: %.4f mm (8阶)\n', Ishii_RMS_Std(end));
fprintf('平均态收敛比率 (8阶/1阶): %.2f%%\n', (Ishii_RMS_Avg(end)/Ishii_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (8阶/1阶): %.2f%%\n', (Ishii_RMS_Std(end)/Ishii_RMS_Std(1))*100);

% 保存统计结果
fprintf('\n保存统计结果...\n');
StatsFile = fullfile('D:\work\Task_Convergence', 'Statistics_Analysis.mat');
save(StatsFile, ...
    'EN4_RMS_Avg', 'EN4_RMS_Std', ...
    'IAP_RMS_Avg', 'IAP_RMS_Std', ...
    'Ishii_RMS_Avg', 'Ishii_RMS_Std', ...
    '-v7.3');
fprintf('统计结果已保存至: %s\n', StatsFile);

% 分析趋势数据
fprintf('\n\n========== 趋势数据统计分析 ==========\n');
trend_dir = 'D:\work\Task_Convergence\Trend_Results';

% 加载趋势数据
fprintf('加载趋势数据...\n');

% EN4 平均态趋势
if exist(fullfile(trend_dir, 'EN4_Average_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'EN4_Average_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    EN4_Trend_TSLA_Avg = trend_TSLA;
    EN4_Trend_SSLA_Avg = trend_SSLA;
    EN4_Trend_HSLA_Avg = trend_HSLA;
    EN4_Trend_Cross_Avg = trend_Cross;
end

% EN4 标准态趋势
if exist(fullfile(trend_dir, 'EN4_StdRef_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'EN4_StdRef_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    EN4_Trend_TSLA_Std = trend_TSLA;
    EN4_Trend_SSLA_Std = trend_SSLA;
    EN4_Trend_HSLA_Std = trend_HSLA;
    EN4_Trend_Cross_Std = trend_Cross;
end

% IAP 平均态趋势
if exist(fullfile(trend_dir, 'IAP_Average_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'IAP_Average_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    IAP_Trend_TSLA_Avg = trend_TSLA;
    IAP_Trend_SSLA_Avg = trend_SSLA;
    IAP_Trend_HSLA_Avg = trend_HSLA;
    IAP_Trend_Cross_Avg = trend_Cross;
end

% IAP 标准态趋势
if exist(fullfile(trend_dir, 'IAP_StdRef_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'IAP_StdRef_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    IAP_Trend_TSLA_Std = trend_TSLA;
    IAP_Trend_SSLA_Std = trend_SSLA;
    IAP_Trend_HSLA_Std = trend_HSLA;
    IAP_Trend_Cross_Std = trend_Cross;
end

% Ishii 平均态趋势
if exist(fullfile(trend_dir, 'Ishii_Average_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'Ishii_Average_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    Ishii_Trend_TSLA_Avg = trend_TSLA;
    Ishii_Trend_SSLA_Avg = trend_SSLA;
    Ishii_Trend_HSLA_Avg = trend_HSLA;
    Ishii_Trend_Cross_Avg = trend_Cross;
end

% Ishii 标准态趋势
if exist(fullfile(trend_dir, 'Ishii_StdRef_Trends.mat'), 'file')
    load(fullfile(trend_dir, 'Ishii_StdRef_Trends.mat'), 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross');
    Ishii_Trend_TSLA_Std = trend_TSLA;
    Ishii_Trend_SSLA_Std = trend_SSLA;
    Ishii_Trend_HSLA_Std = trend_HSLA;
    Ishii_Trend_Cross_Std = trend_Cross;
end

% 分析平均态趋势
fprintf('\n【平均态趋势统计 (单位: mm/年)】\n');
fprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n', ...
    '类别', '最小', '最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值', '标准差');

% EN4 趋势
fprintf('\n--- EN4 趋势 ---.\n');
for order = 1:3
    fprintf('EN4 TSLA趋势 %d阶\t', order);
    stats = calculate_statistics(EN4_Trend_TSLA_Avg(:, :, order));
    data = EN4_Trend_TSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('EN4 SSLA趋势 %d阶\t', order);
    stats = calculate_statistics(EN4_Trend_SSLA_Avg(:, :, order));
    data = EN4_Trend_SSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('EN4 HSLA趋势 %d阶\t', order);
    stats = calculate_statistics(EN4_Trend_HSLA_Avg(:, :, order));
    data = EN4_Trend_HSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
end

% IAP 趋势
fprintf('\n--- IAP 趋势 ---.\n');
for order = 1:3
    fprintf('IAP TSLA趋势 %d阶\t', order);
    stats = calculate_statistics(IAP_Trend_TSLA_Avg(:, :, order));
    data = IAP_Trend_TSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('IAP SSLA趋势 %d阶\t', order);
    stats = calculate_statistics(IAP_Trend_SSLA_Avg(:, :, order));
    data = IAP_Trend_SSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('IAP HSLA趋势 %d阶\t', order);
    stats = calculate_statistics(IAP_Trend_HSLA_Avg(:, :, order));
    data = IAP_Trend_HSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
end

% Ishii 趋势
fprintf('\n--- Ishii 趋势 ---.\n');
for order = 1:3
    fprintf('Ishii TSLA趋势 %d阶\t', order);
    stats = calculate_statistics(Ishii_Trend_TSLA_Avg(:, :, order));
    data = Ishii_Trend_TSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('Ishii SSLA趋势 %d阶\t', order);
    stats = calculate_statistics(Ishii_Trend_SSLA_Avg(:, :, order));
    data = Ishii_Trend_SSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
    
    fprintf('Ishii HSLA趋势 %d阶\t', order);
    stats = calculate_statistics(Ishii_Trend_HSLA_Avg(:, :, order));
    data = Ishii_Trend_HSLA_Avg(:, :, order);
    data = data(~isnan(data));
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean, std(data));
end

% 保存趋势统计结果
fprintf('\n保存趋势统计结果...\n');
TrendStatsFile = fullfile('D:\work\Task_Convergence', 'Trend_Statistics_Analysis.mat');
save(TrendStatsFile, ...
    'EN4_Trend_TSLA_Avg', 'EN4_Trend_SSLA_Avg', 'EN4_Trend_HSLA_Avg', 'EN4_Trend_Cross_Avg', ...
    'EN4_Trend_TSLA_Std', 'EN4_Trend_SSLA_Std', 'EN4_Trend_HSLA_Std', 'EN4_Trend_Cross_Std', ...
    'IAP_Trend_TSLA_Avg', 'IAP_Trend_SSLA_Avg', 'IAP_Trend_HSLA_Avg', 'IAP_Trend_Cross_Avg', ...
    'IAP_Trend_TSLA_Std', 'IAP_Trend_SSLA_Std', 'IAP_Trend_HSLA_Std', 'IAP_Trend_Cross_Std', ...
    'Ishii_Trend_TSLA_Avg', 'Ishii_Trend_SSLA_Avg', 'Ishii_Trend_HSLA_Avg', 'Ishii_Trend_Cross_Avg', ...
    'Ishii_Trend_TSLA_Std', 'Ishii_Trend_SSLA_Std', 'Ishii_Trend_HSLA_Std', 'Ishii_Trend_Cross_Std', ...
    '-v7.3');
fprintf('趋势统计结果已保存至: %s\n', TrendStatsFile);

fprintf('\n>>> 🎉 三套数据统计分析完成! <<<\n');

% 定义统计函数
function stats = calculate_statistics(data)
    % 移除 NaN 值
    valid_data = data(~isnan(data(:)));
    
    % 计算统计量
    stats.min = min(valid_data);
    stats.max = max(valid_data);
    stats.mean = mean(valid_data);
    stats.median = median(valid_data);
    stats.q1 = prctile(valid_data, 25);
    stats.q3 = prctile(valid_data, 75);
    stats.p25 = prctile(valid_data, 2.5);
    stats.p975 = prctile(valid_data, 97.5);
end
