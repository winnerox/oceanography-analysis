%% Test_Matrix_Performance.m
% 功能: 测试 180x360 网格上的随机矩阵计算性能
clear; clc;

% 1. 初始化引擎
fprintf('正在启动引擎...\n');
engine = TEOS10_HighOrder_Engine();

% 2. 生成测试数据 (180x360)
rows = 180;
cols = 360;
fprintf('\n==========================================================\n');
fprintf('   矩阵性能测试 (Matrix Performance Test)\n');
fprintf('   尺寸: %d x %d (共 %d 个点)\n', rows, cols, rows*cols);
fprintf('==========================================================\n');

% 生成"合理"的随机数据
% 盐度 SA: 分布在 30 ~ 40 g/kg 之间
fprintf('生成随机盐度矩阵 (30~40 g/kg)...\n');
SA = 30 + 10 * rand(rows, cols);

% 温度 CT: 分布在 0 ~ 30 度之间
fprintf('生成随机温度矩阵 (0~30 C)...\n');
CT = 0 + 30 * rand(rows, cols);

% 压力 p: 设为标量 1000 dbar (自动扩充)
p = 1000;

% 3. 执行计算 (例如: 4阶混合导数 2温2盐)
n_T = 2; 
n_S = 2;
fprintf('正在计算混合导数 d^4 rho / dCT^2 dSA^2 ...\n');

tic; % 开始计时
val_matrix = engine.calculate_mixed(SA, CT, p, n_T, n_S);
time_cost = toc; % 结束计时

% 4. 结果分析
fprintf('\n----------------------------------------------------------\n');
fprintf('计算耗时 : %.4f 秒\n', time_cost);
fprintf('输出维度 : %d x %d\n', size(val_matrix, 1), size(val_matrix, 2));
fprintf('数据状态 : %s\n', mat2str(size(val_matrix)));

% 检查是否存在 NaN 或 Inf
if any(isnan(val_matrix(:))) || any(isinf(val_matrix(:)))
    fprintf('警告: 结果中包含无效值 (NaN/Inf)！\n');
else
    fprintf('检查通过: 所有数值均有效。\n');
end

% 随机抽查一个点
r_idx = randi(rows);
c_idx = randi(cols);
fprintf('\n随机抽查点 (%d, %d):\n', r_idx, c_idx);
fprintf('  SA = %.4f, CT = %.4f\n', SA(r_idx, c_idx), CT(r_idx, c_idx));
fprintf('  Value = %e\n', val_matrix(r_idx, c_idx));

% 5. 可视化结果 (热力图)
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
imagesc(val_matrix);
colorbar;
title(sprintf('混合导数计算结果 (180x360) - 耗时 %.4fs', time_cost));
xlabel('Grid X'); ylabel('Grid Y');
axis equal tight;