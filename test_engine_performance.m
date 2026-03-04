%% 引擎性能测试脚本
% 功能: 比较 TEOS10_HighOrder_Engine 和 TEOS10_General_Engine 的计算速度
% 测试内容: 引擎初始化、系数计算（1-8阶）
% 测试环境: 随机生成的SA, CT, p数据

clear; clc;

%% 1. 测试配置
max_order = 8;
num_points = 10000; % 测试数据点数
use_cache = false; % 是否使用缓存

%% 2. 生成测试数据
fprintf('>> 生成测试数据...\n');
rng(42); % 设置随机种子，保证结果可重复

% 生成合理范围的测试数据
SA = 30 + 5*rand(num_points, 1);   % 盐度：30-35 psu
CT = 0 + 25*rand(num_points, 1);   % 保守温度：0-25 °C
p = 1000 + 10000*rand(num_points, 1); % 压力：1000-11000 dbar

fprintf('   测试数据：%d个点\n', num_points);
fprintf('   SA范围：%.2f - %.2f psu\n', min(SA), max(SA));
fprintf('   CT范围：%.2f - %.2f °C\n', min(CT), max(CT));
fprintf('   p范围：%.0f - %.0f dbar\n', min(p), max(p));

%% 3. 测试 TEOS10_HighOrder_Engine
fprintf('\n=============================================\n');
fprintf('>> 测试 TEOS10_HighOrder_Engine...\n');
fprintf('=============================================\n');

total_time_ho = tic;

% 初始化引擎
t_ho_init = tic;
if use_cache
    cache_file = 'TEOS10_cache_test.mat';
    engine_ho = TEOS10_HighOrder_Engine(max_order, cache_file);
else
    engine_ho = TEOS10_HighOrder_Engine(max_order);
end
dt_ho_init = toc(t_ho_init);
fprintf('   初始化时间：%.2f 秒\n', dt_ho_init);

% 计算系数
t_ho_calc = tic;
result_ho = engine_ho.calculate_all_orders(SA, CT, p, max_order);
dt_ho_calc = toc(t_ho_calc);
fprintf('   计算时间：%.2f 秒\n', dt_ho_calc);

total_time_ho = toc(total_time_ho);
fprintf('   总时间：%.2f 秒\n', total_time_ho);

%% 4. 测试 TEOS10_General_Engine
fprintf('\n=============================================\n');
fprintf('>> 测试 TEOS10_General_Engine...\n');
fprintf('=============================================\n');

total_time_general = tic;

% 初始化引擎
t_general_init = tic;
engine_general = TEOS10_General_Engine();
dt_general_init = toc(t_general_init);
fprintf('   初始化时间：%.2f 秒\n', dt_general_init);

% 计算系数
t_general_calc = tic;
result_general = engine_general.calculate_all_orders(SA, CT, p, max_order);
dt_general_calc = toc(t_general_calc);
fprintf('   计算时间：%.2f 秒\n', dt_general_calc);

total_time_general = toc(total_time_general);
fprintf('   总时间：%.2f 秒\n', total_time_general);

%% 5. 结果比较
fprintf('\n=============================================\n');
fprintf('>> 性能比较结果\n');
fprintf('=============================================\n');

fprintf('\n=== 时间比较 (秒) ===\n');
fprintf('| 指标 | HighOrder | General | 差异 (%%) |\n');
fprintf('|------|-----------|---------|----------|\n');

% 初始化时间比较
init_diff = ((dt_general_init - dt_ho_init) / dt_ho_init) * 100;
fprintf('| 初始化 | %9.2f | %7.2f | %8.1f%% |\n', dt_ho_init, dt_general_init, init_diff);

% 计算时间比较
calc_diff = ((dt_general_calc - dt_ho_calc) / dt_ho_calc) * 100;
fprintf('| 计算   | %9.2f | %7.2f | %8.1f%% |\n', dt_ho_calc, dt_general_calc, calc_diff);

% 总时间比较
total_diff = ((total_time_general - total_time_ho) / total_time_ho) * 100;
fprintf('| 总计   | %9.2f | %7.2f | %8.1f%% |\n', total_time_ho, total_time_general, total_diff);

% 效率比较
if total_time_ho < total_time_general
    faster_engine = 'TEOS10_HighOrder_Engine';
    speedup = total_time_general / total_time_ho;
    fprintf('\n=== 结论 ===\n');
    fprintf('%s 更快，速度提升 %.2f 倍\n', faster_engine, speedup);
else
    faster_engine = 'TEOS10_General_Engine';
    speedup = total_time_ho / total_time_general;
    fprintf('\n=== 结论 ===\n');
    fprintf('%s 更快，速度提升 %.2f 倍\n', faster_engine, speedup);
end

%% 6. 结果验证
fprintf('\n=== 结果验证 ===\n');

% 比较第一阶导数结果
rel_diff = abs(result_ho.d1_T - result_general.d1_T) ./ mean(abs(result_ho.d1_T));
max_rel_diff = max(rel_diff);
mean_rel_diff = mean(rel_diff);

fprintf('d1_T 相对误差：\n');
fprintf('   最大：%.2e\n', max_rel_diff);
fprintf('   平均：%.2e\n', mean_rel_diff);

if max_rel_diff < 1e-6
    fprintf('✅ 结果一致，相对误差 < 1e-6\n');
else
    fprintf('⚠️  结果存在差异，相对误差 = %.2e\n', max_rel_diff);
end

%% 7. 测试报告
fprintf('\n=============================================\n');
fprintf('>> 测试报告\n');
fprintf('=============================================\n');
fprintf('测试配置：\n');
fprintf('   最大阶数：%d\n', max_order);
fprintf('   测试点数：%d\n', num_points);
fprintf('   是否使用缓存：%s\n', ternary(use_cache, '是', '否'));
fprintf('\n最快引擎：%s\n', faster_engine);
fprintf('速度提升：%.2f 倍\n', speedup);
fprintf('结果一致性：%s\n', ternary(max_rel_diff < 1e-6, '良好', '存在差异'));
fprintf('=============================================\n');

%% 辅助函数
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end