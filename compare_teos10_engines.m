%% Compare_TEOS10_Engines.m
% 这是一个用于对比 TEOS10_General_Engine (工程版) 和 TEOS10_HighOrder_Engine (高阶/符号版)
% 准确性和性能的测试脚本。

clc; clear; close all;

%% 1. 初始化设置
fprintf('=======================================================\n');
fprintf('  TEOS-10 引擎对比测试 (Accuracy & Performance Benchmark)\n');
fprintf('=======================================================\n');

% 尝试添加 GSW 工具箱路径 (如果存在)
gsw_path = 'gsw_matlab_v3_06_16';
if exist(gsw_path, 'dir')
    addpath(genpath(gsw_path));
    fprintf('   [Path] 已添加 GSW 工具箱路径: %s\n', gsw_path);
end

%测试参数
num_points = 1000;      % 测试数据点数量
max_order_test = 4;     % 测试导数阶数 (1-4阶)
cache_file = 'teos10_cache.mat'; % 符号引擎缓存

% 生成随机测试数据 (典型的海洋范围)
rng(42); % 固定种子保证可复现
SA_test = 30 + 10 * rand(num_points, 1);    % 绝对盐度 (30-40 g/kg)
CT_test = 0 + 30 * rand(num_points, 1);     % 保守温度 (0-30 deg C)
p_test  = 0 + 5000 * rand(num_points, 1);   % 压力 (0-5000 dbar)

fprintf('1. 生成测试数据: %d 个点\n', num_points);
fprintf('   SA range: [%.2f, %.2f]\n', min(SA_test), max(SA_test));
fprintf('   CT range: [%.2f, %.2f]\n', min(CT_test), max(CT_test));
fprintf('   p  range: [%.0f, %.0f]\n', min(p_test), max(p_test));
fprintf('\n');

%% 2. 初始化引擎

fprintf('2. 初始化引擎...\n');

% --- 初始化 General Engine (数值版) ---
fprintf('   [General Engine] 初始化... ');
t_gen_init = tic;
try
    gen_engine = TEOS10_General_Engine(max_order_test);
    fprintf('完成 (%.4f s)\n', toc(t_gen_init));
catch ME
    fprintf('失败\n   ❌ %s\n', ME.message);
    return;
end

% --- 初始化 HighOrder Engine (符号版) ---
fprintf('   [HighOrder Engine] 初始化... \n');
% 注意：如果是第一次运行，这里会比较慢，因为它需要符号推导
% 我们开启 use_simplify=false 以加快预编译速度（除非追求极致公式简化）
t_high_init = tic;
try
    high_engine = TEOS10_HighOrder_Engine(max_order_test, cache_file, false);
    fprintf('   [HighOrder Engine] 初始化完成 (%.4f s)\n', toc(t_high_init));
catch ME
    fprintf('   [HighOrder Engine] 初始化失败\n   ❌ %s\n', ME.message);
    fprintf('   ⚠️ 请确保安装了 Symbolic Math Toolbox 且 gsw_specvol.m 可见\n');
    return;
end

fprintf('\n');

%% 3. 运行计算与性能对比

fprintf('3. 运行性能基准测试 (Max Order = %d)...\n', max_order_test);

% --- 运行 General Engine ---
tic;
res_gen = gen_engine.calculate_all_orders(SA_test, CT_test, p_test, max_order_test);
t_gen_run = toc;
fprintf('   [General Engine]   耗时: %.6f s (%.2f μs/pt)\n', t_gen_run, (t_gen_run/num_points)*1e6);

% --- 运行 HighOrder Engine ---
tic;
res_high = high_engine.calculate_all_orders(SA_test, CT_test, p_test, max_order_test);
t_high_run = toc;
fprintf('   [HighOrder Engine] 耗时: %.6f s (%.2f μs/pt)\n', t_high_run, (t_high_run/num_points)*1e6);

ratio = t_high_run / t_gen_run;
fprintf('   >> General Engine 速度是 HighOrder Engine 的 %.1f 倍\n', ratio);
fprintf('\n');

%% 4. 准确性对比 (Accuracy Validation)

fprintf('4. 准确性验证 (以 HighOrder Engine 为基准)...\n');
fprintf('---------------------------------------------------------------------------------\n');
fprintf('| 导数项 |  Max Abs Error  |  Max Rel Error  |          状态           |\n');
fprintf('---------------------------------------------------------------------------------\n');

all_passed = true;
tolerance_rel = 1e-10; % 相对误差容忍度 (数值方法通常能达到机器精度附近)

% 定义要对比的字段
fields = fieldnames(res_gen);

for i = 1:length(fields)
    f_name = fields{i};
    
    if ~isfield(res_high, f_name)
        continue;
    end
    
    val_gen = res_gen.(f_name);
    val_high = res_high.(f_name);
    
    % 计算误差
    abs_err = abs(val_gen - val_high);
    max_abs_err = max(abs_err);
    
    % 计算相对误差 (处理分母接近0的情况)
    % 仅在 val_high 足够大时计算相对误差
    mask = abs(val_high) > 1e-15;
    rel_err = zeros(size(val_gen));
    rel_err(mask) = abs_err(mask) ./ abs(val_high(mask));
    max_rel_err = max(rel_err);
    
    % 判定通过
    if max_rel_err < tolerance_rel
        status_str = '✅ PASS';
    else
        status_str = '⚠️ WARN';
        all_passed = false;
    end
    
    % 这种格式化打印确保列对齐
    fprintf('| %-6s | %15.5e | %15.5e | %s |\n', ...
        f_name, max_abs_err, max_rel_err, status_str);
end

fprintf('---------------------------------------------------------------------------------\n');

if all_passed
    fprintf('\n🎉 结论: 两个引擎结果高度一致！(误差 < 1e-10)\n');
    fprintf('   => General_Engine 可以安全地在生产环境中替代 HighOrder_Engine。\n');
else
    fprintf('\n⚠️ 结论: 发现部分误差较大，请检查实现细节。\n');
end
