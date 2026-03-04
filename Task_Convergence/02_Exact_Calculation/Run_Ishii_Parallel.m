%% Run_Ishii_Parallel.m
% 启动两个独立的 MATLAB 进程并行运行
clear; clc;
disp('启动并行计算任务...');

script1 = 'd:\work\Task_Convergence\02_Exact_Calculation\Calc_Exact_Ishii_Avg.m';
script2 = 'd:\work\Task_Convergence\02_Exact_Calculation\Calc_Exact_Ishii_Std.m';

matlab_exe = 'D:\MATLAB 2023\bin\matlab.exe';

cmd1 = sprintf('"%s" -batch "run(''%s'')" &', matlab_exe, script1);
cmd2 = sprintf('"%s" -batch "run(''%s'')" &', matlab_exe, script2);

fprintf('启动 Ishii_Avg...\n');
system(cmd1, false, 'async');

fprintf('启动 Ishii_Std...\n');
system(cmd2, false, 'async');

disp('两个任务已在后台启动，请等待完成。');
disp('输出文件将保存到: D:\work\Ishii_TSLA_Terms\');
