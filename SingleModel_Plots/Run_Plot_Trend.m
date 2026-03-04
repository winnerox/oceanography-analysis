% 运行整合趋势绘图函数并保持图形窗口打开
% 用法：在MATLAB命令窗口中运行此脚本

clear; clc;

fprintf('========== 运行整合趋势绘图函数 ==========\n');

% 选择要绘制的数据集和状态
% 可以修改这些参数来绘制不同的组合
dataset_name = 'EN4';  % 可选: 'EN4', 'IAP', 'Ishii'
state = 'StdRef';      % 可选: 'StdRef', 'Average'

fprintf('\n绘制 %s (%s) 的趋势图...\n', dataset_name, state);
try
    Plot_Trend_Combined(dataset_name, state);
    fprintf('✓ 绘图成功！图形窗口已打开，请查看。\n');
    fprintf('提示：关闭图形窗口后，MATLAB会继续执行。\n');
    % 等待用户关闭图形窗口
    waitfor(gcf);
catch ME
    fprintf('✗ 绘图失败: %s\n', ME.message);
end

fprintf('\n========== 运行完成 ==========\n');
