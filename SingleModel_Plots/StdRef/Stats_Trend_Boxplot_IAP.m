%% 趋势箱线图 [StdRef Version]
%% Stats_Trend_Boxplot_IAP.m
% 功能: IAP 数据趋势箱线图 (8 阶展开项) - 单轴修正版
clear; clc; close all;

%% 0. 配置
OutputDir = 'D:\work\Figures';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 1. 加载数据
files = {'IAP_TSLA_Terms_1to8_StdRef.mat'};
input_file = '';
DataDir = 'D:\work\MAT_Data';
for i = 1:length(files)
    if exist(files{i}, 'file')
        input_file = files{i};
    elseif exist(fullfile(DataDir, files{i}), 'file')
        input_file = fullfile(DataDir, files{i});
    end
    if ~isempty(input_file), break; end
end
if isempty(input_file), error('❌ 找不到 IAP 数据文件: %s', files{1}); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 变量名统一处理
if ~exist('time_axis', 'var') && exist('Time_Axis', 'var')
    time_axis = Time_Axis;
end

% 准备趋势数据
if exist('TSLA_AllOrders','var')
    slice1 = TSLA_AllOrders(:,:,1,1);
elseif exist('TSLA_Result','var')
    TSLA_AllOrders = TSLA_Result;
end

TSLA_AllOrders = TSLA_AllOrders(:, :, :, :);
mask = ~isnan(TSLA_AllOrders(:,:,1,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Trend = nan(num_valid, 8);
t_vec = time_axis(:);

for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    vec_t = nan(num_valid, 1);
    for k = 1:num_valid
        ts = squeeze(slice(rows(k), cols(k), :));
        if any(isnan(ts)), continue; end
        p = polyfit(t_vec, ts, 1);
        vec_t(k) = p(1); 
    end
    Data_Trend(:,n) = vec_t;
end

% 智能检测单位
if max(abs(Data_Trend(:))) < 0.5 && max(abs(Data_Trend(:))) > 0
    fprintf('>> 检测到数值过小，推测为 m/yr，自动 * 1000\n');
    Data_Trend = Data_Trend * 1000;
end

%% 绘图参数
BoxBgColor = [0.95 0.95 1];
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

% 颜色定义
C1 = [0.8500 0.3250 0.0980]; % Red
C2 = [0 0.4470 0.7410];      % Blue
C3 = [0.4660 0.6740 0.1880]; % Green

% 箱线图位置
%% 绘图 - 0轴对齐版
fprintf('>> 绘制 IAP 趋势箱线图 (0轴对齐版)...\n');
figure('Position', [50, 50, 1200, 600], 'Color', 'w', 'Name', 'IAP Trend ZeroAligned');

% --- 绘图参数 ---
MainPos = [0.07, 0.14, 0.80, 0.80];
BoxWidth = 0.6;

% 定义每个组的中心位置
Pos_G1 = [1, 2, 3];
Pos_G2 = [4.2, 5.2];
Pos_G3 = [6.4, 7.4, 8.4];

X_Limit = [0.4, 9.0];

% --- Layer 1 (Red) ---
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:3
    draw_jitter_unit(ax1, Pos_G1(k), Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
% 参考振幅图逻辑: 
max_v1 = max(abs(prctile(Data_Trend(:,1:3), [0.2 99.8], 'all'))); 
if isnan(max_v1)||max_v1==0, max_v1=1; end
Limit1 = max_v1 * 1.1; 
% 关键: 只需要设置 ylim 对称即可，不需要强制 tick
ylim(ax1, [-Limit1, Limit1]); 

set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% --- Layer 2 (Blue) ---
ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 1:2
    draw_jitter_unit(ax2, Pos_G2(k), Data_Trend(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
end
max_v2 = max(abs(prctile(Data_Trend(:,4:5), [0.2 99.8], 'all'))); 
if isnan(max_v2)||max_v2==0, max_v2=1; end
Limit2 = max_v2 * 1.2;
ylim(ax2, [-Limit2, Limit2]); 

set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% --- Layer 3 (Green) ---
Offset_Ratio = 0.08; 
Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
ax3 = axes('Position', Pos3); hold(ax3, 'on');
Scale_Factor = Pos3(3) / MainPos(3);
Real_X_Limit = X_Limit * Scale_Factor; 

for k = 1:3
    draw_jitter_unit(ax3, Pos_G3(k), Data_Trend(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
end
max_v3 = max(abs(prctile(Data_Trend(:,6:8), [0.2 99.8], 'all'))); 
if isnan(max_v3)||max_v3==0, max_v3=1; end
Limit3 = max_v3 * 1.2;
ylim(ax3, [-Limit3, Limit3]);

set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% --- 装饰与分割线 ---
yline(ax1, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.2); 
xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels);
xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, '\fontname{Microsoft YaHei}IAP (StdRef) 各阶趋势分布图 (0轴对齐版)', 'FontSize', 16, 'FontWeight', 'bold');

% === 【图例修正】Figure 1 (右上角) ===
h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out  = scatter(ax1, NaN, NaN, ScatterSize, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.6);


legend(ax1, [h_box, h_med, h_mean, h_out], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', ...      % 右上角
       'FontSize', 10, ...
       'Box', 'on', ...                  % 保留白色背景
       'Color', [1 1 1], ...             % 纯白
       'EdgeColor', 'none');             % 无边框

% === 交互优化：锁定Y轴缩放，同步X轴 ===
% === 交互优化：多轴同步缩放 (X & Y) ===
Ratio2 = Limit2 / Limit1;
Ratio3 = Limit3 / Limit1;

linkaxes([ax1, ax2], 'x'); 
hZoom = zoom(gcf);
hZoom.Motion = 'both'; 
hZoom.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
hZoom.Enable = 'on';

hPan = pan(gcf);
hPan.Motion = 'both';
hPan.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);






%% ==========================================================
%% Figure 2: 单 Y 轴 (紧凑版)
%% ==========================================================
fprintf('>> 绘制 Figure 2: 紧凑版 单 Y 轴趋势图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Trend Single-Compact');

% 【关键】创建坐标轴 (与 EN4 保持一致)
ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');

% 统一颜色 (IAP Blue)
UniColor = [0 0.4470 0.7410]; 

for i = 1:8
    draw_jitter_unit(ax, i, Data_Trend(:,i), 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha);
end

yline(ax, 0, '--', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
set(ax, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman', 'XColor', 'k', 'YColor', 'k');
All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
set(ax, 'XTick', 1:8, 'XTickLabel', All_Labels);
xlabel(ax, '\fontname{Microsoft YaHei}泰勒展开项', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax, '\fontname{Microsoft YaHei}趋势变化率 (mm/yr)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-10, 10]); 
grid on; ax.GridAlpha = 0.15;
title('\fontname{Microsoft YaHei}IAP (StdRef) 各阶趋势分布图 (单轴)', 'FontSize', 16, 'FontWeight', 'bold');

% === 【图例修正】Figure 2 (右上角) ===
h_box2  = patch(ax, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med2  = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean2 = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out2  = scatter(ax, NaN, NaN, ScatterSize, UniColor, 'filled', 'MarkerFaceAlpha', 0.6);

legend(ax, [h_box2, h_med2, h_mean2, h_out2], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', ...      % 右上角
       'FontSize', 10, ...
       'Box', 'on', ...                  % 保留白色背景
       'Color', [1 1 1], ...             % 纯白
       'EdgeColor', 'none');             % 无边框

fprintf('>> 绘图完成！\n');

%% 辅助函数
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w = prctile(data, 97.5);
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    
    if ~isempty(outliers)
        if length(outliers) > 2000
            rand_idx = randperm(length(outliers), 2000);
            outliers = outliers(rand_idx);
        end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        try
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s);
        end
    end
    
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

% 辅助函数: 全轴同步 (X & Y) + 强制 0 轴对齐
function SyncAllAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        % 1. 强制主轴 Y 范围对称 (0 轴居中)
        L = max(abs(ax1.YLim));
        ax1.YLim = [-L, L];
        
        % 2. 依次同步次轴 (保持比例)
        ax2.YLim = [-L*R2, L*R2];
        ax3.YLim = [-L*R3, L*R3];

        % 3. 同步 ax3 的 X 轴
        ax3.XLim = ax1.XLim * Factor;
    catch
    end
end
