%% Analyze_Statistics_Formatted.m
% =========================================================================
% 功能：分析三套数据集（EN4、IAP、Ishii）的泰勒展开收敛性统计信息
% 输出：格式化的统计分析结果，包括美观的表格和汇总
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

% 加载 EN4 数据
load(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
load(fullfile(EN4_DataDir, 'EN4_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
EN4_TSLA = TSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
EN4_SSLA = SSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
EN4_HSLA = EN4_TSLA + EN4_SSLA + Cross_AllOrders;

% 加载 IAP 数据
load(fullfile(IAP_DataDir, 'IAP_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
IAP_Exact_Avg = TSLA_Exact_Avg;
load(fullfile(IAP_DataDir, 'IAP_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
IAP_TSLA = TSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
IAP_SSLA = SSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
IAP_HSLA = IAP_TSLA + IAP_SSLA + Cross_AllOrders;

% 加载 Ishii 数据
load(fullfile(Ishii_DataDir, 'Ishii_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg');
Ishii_Exact_Avg = TSLA_Exact_Avg;
load(fullfile(Ishii_DataDir, 'Ishii_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
Ishii_TSLA = TSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders');
Ishii_SSLA = SSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders');
Ishii_HSLA = Ishii_TSLA + Ishii_SSLA + Cross_AllOrders;

% 主输出
fprintf('\n=================================================================\n');
fprintf('                    三套数据泰勒展开统计分析\n');
fprintf('=================================================================\n\n');

% 分析低阶项（1-3阶）
fprintf('【子图1：第1-3阶 - SSLA, TSLA, HSLA】\n\n');

% EN4 低阶
fprintf('--- EN4 低阶 (1-3阶) ---.\n\n');
print_table_header();
for order = 1:3
    stats = calculate_statistics(EN4_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 3
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% IAP 低阶
fprintf('\n--- IAP 低阶 (1-3阶) ---.\n\n');
print_table_header();
for order = 1:3
    stats = calculate_statistics(IAP_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 3
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% Ishii 低阶
fprintf('\n--- Ishii 低阶 (1-3阶) ---.\n\n');
print_table_header();
for order = 1:3
    stats = calculate_statistics(Ishii_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 3
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% 分析中阶项（4-6阶）
fprintf('\n\n【子图2：第4-6阶 - SSLA, TSLA, HSLA】\n\n');

% EN4 中阶
fprintf('--- EN4 中阶 (4-6阶) ---.\n\n');
print_table_header();
for order = 4:6
    stats = calculate_statistics(EN4_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 6
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% IAP 中阶
fprintf('\n--- IAP 中阶 (4-6阶) ---.\n\n');
print_table_header();
for order = 4:6
    stats = calculate_statistics(IAP_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 6
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% Ishii 中阶
fprintf('\n--- Ishii 中阶 (4-6阶) ---.\n\n');
print_table_header();
for order = 4:6
    stats = calculate_statistics(Ishii_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 6
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% 分析高阶项（7-8阶）
fprintf('\n\n【子图3：第7-8阶 - SSLA, TSLA, HSLA】\n\n');

% EN4 高阶
fprintf('--- EN4 高阶 (7-8阶) ---.\n\n');
print_table_header();
for order = 7:8
    stats = calculate_statistics(EN4_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(EN4_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 8
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% IAP 高阶
fprintf('\n--- IAP 高阶 (7-8阶) ---.\n\n');
print_table_header();
for order = 7:8
    stats = calculate_statistics(IAP_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(IAP_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 8
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% Ishii 高阶
fprintf('\n--- Ishii 高阶 (7-8阶) ---.\n\n');
print_table_header();
for order = 7:8
    stats = calculate_statistics(Ishii_SSLA(:, :, :, order));
    print_table_row(sprintf('SSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_TSLA(:, :, :, order));
    print_table_row(sprintf('TSLA %d阶', order), stats);
    stats = calculate_statistics(Ishii_HSLA(:, :, :, order));
    print_table_row(sprintf('HSLA %d阶', order), stats);
    if order < 8
        fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    end
end
print_table_footer();

% 收敛性统计汇总
fprintf('\n\n=================================================================\n');
fprintf('                       收敛性统计汇总\n');
fprintf('=================================================================\n\n');

% 打印收敛性表格
fprintf('+----------+----------------+----------------+------------------+------------------+\n');
fprintf('| 数据集   | 平均态最终残差 | 标准态最终残差 | 平均态收敛比率   | 标准态收敛比率   |\n');
fprintf('|          |     (mm)       |     (mm)       |     (8阶/1阶)    |     (8阶/1阶)    |\n');
fprintf('+----------+----------------+----------------+------------------+------------------+\n');

% EN4 收敛性
fprintf('| EN4      | %14.4f | %14.4f | %16.2f%% | %16.2f%% |\n', ...
    EN4_RMS_Avg(end), EN4_RMS_Std(end), ...
    (EN4_RMS_Avg(end)/EN4_RMS_Avg(1))*100, (EN4_RMS_Std(end)/EN4_RMS_Std(1))*100);

% IAP 收敛性
fprintf('| IAP      | %14.4f | %14.4f | %16.2f%% | %16.2f%% |\n', ...
    IAP_RMS_Avg(end), IAP_RMS_Std(end), ...
    (IAP_RMS_Avg(end)/IAP_RMS_Avg(1))*100, (IAP_RMS_Std(end)/IAP_RMS_Std(1))*100);

% Ishii 收敛性
fprintf('| Ishii    | %14.4f | %14.4f | %16.2f%% | %16.2f%% |\n', ...
    Ishii_RMS_Avg(end), Ishii_RMS_Std(end), ...
    (Ishii_RMS_Avg(end)/Ishii_RMS_Avg(1))*100, (Ishii_RMS_Std(end)/Ishii_RMS_Std(1))*100);

fprintf('+----------+----------------+----------------+------------------+------------------+\n');

% 保存统计结果
fprintf('\n保存统计结果...\n');
StatsFile = fullfile('D:\work\Task_Convergence', 'Statistics_Analysis_Formatted.mat');
save(StatsFile, ...
    'EN4_RMS_Avg', 'EN4_RMS_Std', ...
    'IAP_RMS_Avg', 'IAP_RMS_Std', ...
    'Ishii_RMS_Avg', 'Ishii_RMS_Std', ...
    '-v7.3');
fprintf('统计结果已保存至: %s\n', StatsFile);

fprintf('\n=================================================================\n');
fprintf('                        分析完成\n');
fprintf('=================================================================\n');

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

% 定义格式化输出函数
function print_table_header()
    fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
    fprintf('| 类别     | 原始最小  | 原始最大  | 剔除后最小 | 剔除后最大 |  2.5%%   |   25%%   |  中位数  |   75%%   |  97.5%%  |  均值    |\n');
    fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
end

function print_table_row(category, stats)
    fprintf('| %-8s | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f |\n', ...
        category, stats.min, stats.max, stats.min, stats.max, stats.p25, stats.q1, stats.median, stats.q3, stats.p975, stats.mean);
end

function print_table_footer()
    fprintf('+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n');
end
