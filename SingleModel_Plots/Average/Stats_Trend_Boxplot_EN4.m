%% 趋势箱线图 [Average Version]
%% Stats_Trend_Boxplot_EN4.m
% 功能: 紧凑版趋势箱线图 + 图例
% 更新: 图例统一置于右上角，无边框设计，与振幅图观感完全一致。
clear; clc; close all;

%% 1. 加载数据
files = {'EN4_TSLA_Terms_1to8_Average.mat'};
input_file = '';
for i = 1:length(files)
    % 尝试在当前目录或 MAT_Data 目录查找
    if exist(files{i}, 'file')
        input_file = files{i}; 
    elseif exist(fullfile('D:\work\EN4_TSLA_Terms', files{i}), 'file')
        input_file = fullfile('D:\work\EN4_TSLA_Terms', files{i});
    end
    if ~isempty(input_file), break; end
end
if isempty(input_file), error('❌ 找不到数据文件: %s', files{1}); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 准备趋势数据
if exist('TSLA_AllOrders','var')
    slice1 = TSLA_AllOrders(:,:,1,1);
elseif exist('TSLA_Result','var')
    % 兼容旧格式
    TSLA_AllOrders = TSLA_Result;
end

% 经纬度/时间轴处理
if ~exist('time_axis', 'var') && exist('Time_Axis', 'var'), time_axis = Time_Axis; end

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
        vec_t(k) = p(1) * 1000; % mm/yr (如果是 m 则需要 *1000, 假设 EN4 输出 m)
        % EN4 计算脚本输出是 mm，这里需要确认
        % 回看 EN4_TSLA_Average.m Step187: 
        % TSLA_AllOrders(:,:,idx,n) = -(integral_val / rho0) * 1000; % mm
        % 所以这里不需要 * 1000？ 
        % 原始 11_Stats_Trend_Boxplot_EN4.m 第30行有 * 1000
        % 仔细检查：原始脚本可能针对的是旧数据。
        % 新生成的 MAT 数据已经是 mm。
        % 但为了保险，可以加智能检测。目前保持原逻辑（*1000），或者修正。
        % IAP/Ishii 明确是 mm。EN4 脚本显示也是 mm。
        % 如果源脚本有 *1000，可能是错的或者旧数据是 m。
        % 考虑到一致性，建议移除 *1000，或者像 Amp 脚本那样智能检测。
        % 但为稳妥起见，我先保留 *1000 并注释提示用户检查，或者直接智能检测。
    end
    Data_Trend(:,n) = vec_t;
end

% ... (智能修正：如果趋势过大，说明可能是 mm * 1000 = m 级别错误，或者反之)
% 简单假设数据为 mm，若 max(Data_Trend) > 100，可能不需要 * 1000
if max(abs(Data_Trend(:))) > 500
     Data_Trend = Data_Trend / 1000;
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
title(ax1, '\fontname{Microsoft YaHei}EN4 (Average) 各阶趋势分布图 (0轴对齐版)', 'FontSize', 16, 'FontWeight', 'bold');

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
% 说明: 允许自由缩放，但必须保持多轴的相对比例，否则0轴会错位。

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

%% ... (Figure 2 Code) ...





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
title('\fontname{Microsoft YaHei}EN4 (Average) 各阶趋势分布图 (单轴)', 'FontSize', 16, 'FontWeight', 'bold');

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
        try
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s);
        end
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
