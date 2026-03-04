%% 振幅箱线图 [Average Version]
%% Stats_Amp_Boxplot_IAP.m
% 功能: IAP 数据振幅箱线图 (8 阶展开项) - 增强版
% 特性: 智能单位检测
clear; clc; close all;

%% 0. 配置
OutputDir = 'D:\work\Figures';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 1. 加载数据
files = {'IAP_TSLA_Terms_1to8_Average.mat'};
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

if ~exist('time_axis', 'var') && exist('Time_Axis', 'var')
    time_axis = Time_Axis;
end

% 准备数据
if exist('TSLA_AllOrders','var')
    slice1 = TSLA_AllOrders(:,:,1,1);
elseif exist('TSLA_Result','var')
    TSLA_AllOrders = TSLA_Result;
end

TSLA_AllOrders = TSLA_AllOrders(:, :, :, :);
mask = ~isnan(TSLA_AllOrders(:,:,1,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Amp = nan(num_valid, 8); 
t_vec = time_axis(:);

%% 2. 智能单位检测与计算
fprintf('>> [智能检测] 检查数据单位量级...\n');
SampleT1 = abs(TSLA_AllOrders(:,:,:,1));
MeanT1 = nanmean(SampleT1(:));
if MeanT1 < 5
    fprintf('   检测到 T1 均值 = %.4f -> 推测单位为 [米]\n', MeanT1);
    UnitFactor = 1000; % 米 -> 毫米
else
    fprintf('   检测到 T1 均值 = %.4f -> 推测单位为 [毫米]\n', MeanT1);
    UnitFactor = 1;    % 已经是毫米
end

fprintf('>> 计算振幅数据 (Residual STD, mm)...\n');
for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    vec_a = nan(num_valid, 1);
    for k = 1:num_valid
        ts = squeeze(slice(rows(k), cols(k), :));
        idx_valid = ~isnan(ts);
        if sum(idx_valid) < length(ts) * 0.5, continue; end
        p = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
        detrended_ts = ts(idx_valid) - polyval(p, t_vec(idx_valid));
        vec_a(k) = std(detrended_ts) * UnitFactor; 
    end
    Data_Amp(:,n) = vec_a;
end

%% 3. 绘图参数
C1 = [0.85 0.325 0.098]; % 红
C2 = [0 0.447 0.741];    % 蓝
C3 = [0.466 0.674 0.188];% 绿 (Ishii 特色可能偏绿，但这里用统一多轴配色)
BoxBgColor = [0.95 0.95 1];
UniColor   = [0 0.447 0.741]; % IAP用统一蓝色 (单轴)

BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.5; ScatterSize = 6; % 减小点大

Pos_G1 = [1, 2];       % T1, T2
Pos_G2 = [3.8, 4.8, 5.8];  % T3, T4, T5 (Shifted slightly for tight layout)
Pos_G3 = [7.5, 8.5, 9.5];  % T6, T7, T8
X_Limit = [0.4, 10.2]; 

%% ==========================================================
%% Figure 1: 3-Y 轴 (紧凑版)
%% ==========================================================
fprintf('>> 绘制 Figure 1: IAP 多轴紧凑版...\n');
f1 = figure('Color', 'w', 'Name', 'IAP Amp Compact', 'WindowState', 'maximized');
MainPos = [0.07, 0.14, 0.80, 0.80]; 

% --- Layer 1 (Red): Dominant (T1, T2) ---
ax1 = axes('Position', MainPos); hold(ax1, 'on'); 
for k = 1:length(Pos_G1)
    draw_jitter_unit(ax1, Pos_G1(k), Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = prctile(Data_Amp(:,1:2), 99.8, 'all'); if isnan(max_v1)||max_v1==0, max_v1=1; end
ylim(ax1, [0, max_v1*1.1]); 
set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% --- Layer 2 (Blue): Secondary (T3, T4, T5) ---
ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 1:length(Pos_G2)
    draw_jitter_unit(ax2, Pos_G2(k), Data_Amp(:,k+2), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
end
max_v2 = prctile(Data_Amp(:,3:5), 99.8, 'all'); if isnan(max_v2)||max_v2==0, max_v2=1; end
ylim(ax2, [0, max_v2*1.2]);
set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% --- Layer 3 (Green): High-Order (T6, T7, T8) ---
Offset_Ratio = 0.08; 
Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
ax3 = axes('Position', Pos3); hold(ax3, 'on');
Scale_Factor = Pos3(3) / MainPos(3);
Real_X_Limit = X_Limit * Scale_Factor; 

for k = 1:length(Pos_G3)
    draw_jitter_unit(ax3, Pos_G3(k), Data_Amp(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
end
max_v3 = prctile(Data_Amp(:,6:8), 99.8, 'all'); if isnan(max_v3)||max_v3==0, max_v3=1; end
ylim(ax3, [0, max_v3*1.2]);
set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% --- 装饰 ---
yline(ax1, 0, '-', 'Color', 'k', 'LineWidth', 1.2); 
xline(ax1, 2.9, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5); % 分割线 1
xline(ax1, 6.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5); % 分割线 2

All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels);
xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
title(ax1, '\fontname{Microsoft YaHei}IAP (Average) 各阶独立项振幅 (多轴布局)', 'FontSize', 18, 'FontWeight', 'bold');

% 图例
h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
h_out  = scatter(ax1, NaN, NaN, ScatterSize, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.6);
legend(ax1, [h_box, h_med, h_out], {'IQR Box (25%-75%)', 'Median', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');

%% ==========================================================
%% Figure 2: 单 Y 轴 (标准版)
%% ==========================================================
fprintf('>> 绘制 Figure 2: IAP 单轴版...\n');
f2 = figure('Color', 'w', 'Name', 'IAP Amp Single', 'WindowState', 'maximized');
axS = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(axS, 'on');

for i = 1:8
    draw_jitter_unit(axS, i, Data_Amp(:,i), 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha);
end

yline(axS, 0, '-', 'Color', 'k', 'LineWidth', 1.5);
set(axS, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');
set(axS, 'XTick', 1:8, 'XTickLabel', All_Labels);
xlabel(axS, '\fontname{Microsoft YaHei}泰勒展开项', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(axS, '\fontname{Microsoft YaHei}振幅强度 (mm)', 'FontSize', 16, 'FontWeight', 'bold');

CurrentMax = max(prctile(Data_Amp, 99.5, 'all'), 1);
ylim([0, CurrentMax * 1.2]); 
grid on; axS.GridAlpha = 0.15;
title('\fontname{Microsoft YaHei}IAP (Average) 各阶独立项振幅 (单轴)', 'FontSize', 18, 'FontWeight', 'bold');

h_box_s  = patch(axS, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med_s  = plot(axS, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean_s = plot(axS, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out_s  = scatter(axS, NaN, NaN, ScatterSize, UniColor, 'filled', 'MarkerFaceAlpha', 0.6);
legend(axS, [h_box_s, h_med_s, h_mean_s, h_out_s], ...
       {'IQR Box', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');

fprintf('>> 全部完成。\n');

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
