%% Step31_Calc_Nonlinearity_Index.m
% 功能:
%   计算温度驱动 TSLA 的非线性强度指标 (NLI)
%   包括:
%     1) 绝对非线性强度 NLI_abs
%     2) 相对非线性强度 NLI_rel (最重要)
%     3) 带符号非线性指数 NLI_signed
%
clear; clc;

%% ================= 用户配置 =================
input_file  = 'EN4_TSLA_Terms_1to15_StdRef.mat';
output_file = 'EN4_TSLA_Nonlinearity_Index.mat';
% ============================================

fprintf('>> 读取高阶 TSLA 数据...\n');
load(input_file, 'TSLA_AllOrders', 'time_axis', 'lat', 'lon');

% TSLA 分解
TSLA_1  = TSLA_AllOrders(:,:,:,1);        % 线性项
TSLA_NL = TSLA_AllOrders(:,:,:,2:end);    % 非线性项

%% 1. 绝对非线性强度 (m)
fprintf('>> 计算绝对非线性强度 NLI_abs...\n');
NLI_abs = nansum(abs(TSLA_NL), 4);

%% 2. 相对非线性强度 (无量纲, 核心指标)
fprintf('>> 计算相对非线性强度 NLI_rel...\n');
NLI_rel = NLI_abs ./ abs(TSLA_1);

%% 3. 带符号非线性指数
fprintf('>> 计算带符号非线性指数 NLI_signed...\n');
NLI_signed = nansum(TSLA_NL, 4) ./ TSLA_1;

%% 4. 清理异常值
NLI_rel(NLI_rel > 5) = NaN;        % 防止极端值污染统计
NLI_signed(abs(NLI_signed) > 5) = NaN;

%% 5. 保存
fprintf('>> 保存非线性指标...\n');
save(output_file, ...
     'NLI_abs', 'NLI_rel', 'NLI_signed', ...
     'time_axis', 'lat', 'lon', '-v7.3');

fprintf('>> 非线性强度指标计算完成!\n');
