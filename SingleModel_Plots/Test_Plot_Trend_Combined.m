% 测试整合趋势绘图函数
% 运行示例：绘制EN4标准态和平均态的趋势图

clear; clc;

fprintf('========== 测试整合趋势绘图函数 ==========\n');

% 测试EN4标准态
fprintf('\n1. 测试 EN4 StdRef...\n');
try
    Plot_Trend_Combined('EN4', 'StdRef');
    fprintf('✓ EN4 StdRef 绘图成功！\n');
catch ME
    fprintf('✗ EN4 StdRef 绘图失败: %s\n', ME.message);
end

% 测试EN4平均态
fprintf('\n2. 测试 EN4 Average...\n');
try
    Plot_Trend_Combined('EN4', 'Average');
    fprintf('✓ EN4 Average 绘图成功！\n');
catch ME
    fprintf('✗ EN4 Average 绘图失败: %s\n', ME.message);
end

% 测试IAP标准态
fprintf('\n3. 测试 IAP StdRef...\n');
try
    Plot_Trend_Combined('IAP', 'StdRef');
    fprintf('✓ IAP StdRef 绘图成功！\n');
catch ME
    fprintf('✗ IAP StdRef 绘图失败: %s\n', ME.message);
end

% 测试IAP平均态
fprintf('\n4. 测试 IAP Average...\n');
try
    Plot_Trend_Combined('IAP', 'Average');
    fprintf('✓ IAP Average 绘图成功！\n');
catch ME
    fprintf('✗ IAP Average 绘图失败: %s\n', ME.message);
end

% 测试Ishii标准态
fprintf('\n5. 测试 Ishii StdRef...\n');
try
    Plot_Trend_Combined('Ishii', 'StdRef');
    fprintf('✓ Ishii StdRef 绘图成功！\n');
catch ME
    fprintf('✗ Ishii StdRef 绘图失败: %s\n', ME.message);
end

% 测试Ishii平均态
fprintf('\n6. 测试 Ishii Average...\n');
try
    Plot_Trend_Combined('Ishii', 'Average');
    fprintf('✓ Ishii Average 绘图成功！\n');
catch ME
    fprintf('✗ Ishii Average 绘图失败: %s\n', ME.message);
end

fprintf('\n========== 测试完成 ==========\n');
