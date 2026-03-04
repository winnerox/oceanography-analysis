%% 趋势箱线图
%% Step47_Trend_Boxplot_Compact_Final.m
% 功能: 紧凑版趋势箱线图 + 图例
% 更新: 图例统一置于右上角，无边框设计，与振幅图观感完全一致。
clear; clc; close all;

%% 1. 加载数据
files = {'EN4_TSLA_Terms_1to15_Safe.mat', 'EN4_TSLA_Terms_1to15_GPU.mat', 'EN4_TSLA_Terms_1to15_StdRef.mat'};
input_file = '';
for i = 1:length(files), if exist(files{i}, 'file'), input_file = files{i}; break; end; end
if isempty(input_file), error('❌ 找不到数据文件'); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 准备趋势数据
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
        vec_t(k) = p(1) * 1000; 
    end
    Data_Trend(:,n) = vec_t;
end

%% ==========================================================
%% 绘图参数
%% ==========================================================
C1 = [0.85 0.325 0.098]; % 红
C2 = [0 0.447 0.741];    % 蓝
C3 = [0.466 0.674 0.188];% 绿
UniColor = [0.2 0.4 0.7]; 
BoxBgColor = [0.95 0.95 1];

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

Pos_G1 = [1, 2, 3];
Pos_G2 = [4.2, 5.2]; 
Pos_G3 = [6.4, 7.4, 8.4];
X_Limit = [0.4, 9.0];

%% ==========================================================
%% Figure 1: 3-Y 轴 (紧凑版)
%% ==========================================================
fprintf('>> 绘制 Figure 1: 紧凑版 3-Y 轴趋势图...\n');
figure('Position', [50, 500, 1200, 550], 'Color', 'w', 'Name', 'Trend Compact');
MainPos = [0.07, 0.14, 0.80, 0.80]; 

% --- Layer 1 (Red) ---
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:3
    draw_jitter_unit(ax1, Pos_G1(k), Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = prctile(abs(Data_Trend(:,1:3)), 99.8, 'all'); 
ylim(ax1, [-max_v1*1.1, max_v1*1.1]);
set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% --- Layer 2 (Blue) ---
ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 1:2
    draw_jitter_unit(ax2, Pos_G2(k), Data_Trend(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
end
max_v2 = prctile(abs(Data_Trend(:,4:5)), 99.8, 'all'); 
ylim(ax2, [-max_v2*1.2, max_v2*1.2]);
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
max_v3 = prctile(abs(Data_Trend(:,6:8)), 99.8, 'all'); 
ylim(ax3, [-max_v3*1.2, max_v3*1.2]);
set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% --- 装饰与分割线 ---
yline(ax1, 0, '--', 'Color', [.6 .6 .6], 'LineWidth', 1.2);
xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels);
xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, '\fontname{Microsoft YaHei}各阶独立项趋势多尺度分布图 (紧凑布局)', 'FontSize', 16, 'FontWeight', 'bold');

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

%% ==========================================================
%% Figure 2: 单 Y 轴 (紧凑版)
%% ==========================================================
fprintf('>> 绘制 Figure 2: 紧凑版 单 Y 轴趋势图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Trend Single-Compact');

% 【关键】创建坐标轴
ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');

for i = 1:8
    draw_jitter_unit(ax, i, Data_Trend(:,i), 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha);
end

yline(ax, 0, '--', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
set(ax, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman', 'XColor', 'k', 'YColor', 'k');
set(ax, 'XTick', 1:8, 'XTickLabel', All_Labels);
xlabel(ax, '\fontname{Microsoft YaHei}泰勒展开项', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax, '\fontname{Microsoft YaHei}趋势变化率 (mm/yr)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-10, 10]); 
grid on; ax.GridAlpha = 0.15;
title('\fontname{Microsoft YaHei}各阶独立项趋势多尺度分布图 (单轴紧凑)', 'FontSize', 16, 'FontWeight', 'bold');

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

%% ============================================================
%% 内部函数
%% ============================================================
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25);
    q3 = prctile(data, 75);
    med_val = median(data);
    mean_val = mean(data);
    low_w = prctile(data, 2.5); 
    up_w  = prctile(data, 97.5);
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    
    if ~isempty(outliers)
        if length(outliers) > 2000
            rand_idx = randperm(length(outliers), 2000);
            outliers = outliers(rand_idx);
        end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
    end
    
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    x_L = center - w_box/2;
    x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, ...
        'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end