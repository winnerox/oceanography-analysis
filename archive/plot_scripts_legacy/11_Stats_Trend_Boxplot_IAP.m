%% Step18_IAP_Trend_Boxplot.m
% =========================================================================
% 功能: IAP 数据趋势箱线图 (8 阶展开项) - 单轴修正版
% =========================================================================
clear; clc; close all;

%% 0. 配置
OutputDir = 'D:\work\Figures';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 1. 加载数据
DataDir = 'D:\work\MAT_Data';
files = {'IAP_TSLA_Terms_1to8_Average.mat', 'IAP_TSLA_Result.mat'};
input_file = '';
for i = 1:length(files)
    fp = fullfile(DataDir, files{i});
    if exist(fp, 'file'), input_file = fp; break; end
end
if isempty(input_file), error('❌ 找不到 IAP 数据文件'); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 变量名统一处理
if ~exist('time_axis', 'var') && exist('Time_Axis', 'var')
    time_axis = Time_Axis;
end

% 准备趋势数据
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
        vec_t(k) = p(1); % 【修正】源数据已经是 mm，不需要再 * 1000
    end
    Data_Trend(:,n) = vec_t;
end

%% 绘图参数
UniColor = [0 0.447 0.741];  % IAP 蓝色
BoxBgColor = [0.95 0.95 1];
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

%% 绘图 - 单轴设计 (保持与 EN4 一致)
fprintf('>> 绘制 IAP 趋势箱线图 (单轴修正)...\n');
figure('Position', [50, 50, 1200, 600], 'Color', 'w', 'Name', 'IAP Trend SingleAxis');
ax = axes('Position', [0.08, 0.12, 0.90, 0.80]); hold(ax, 'on');

All_Labels = {'T1','T2','T3','T4','T5','T6','T7','T8'};

for i = 1:8
    draw_jitter_unit(ax, i, Data_Trend(:,i), 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha);
end

yline(ax, 0, '--', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

set(ax, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
set(ax, 'XTick', 1:8, 'XTickLabel', All_Labels);
xlabel(ax, '\fontname{Microsoft YaHei}泰勒展开项', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax, '\fontname{Microsoft YaHei}趋势变化率 (mm/yr)', 'FontSize', 14, 'FontWeight', 'bold');

% 【关键】为了与 EN4 观感一致，放开 Y 轴或设置类似的范围
% EN4 范围大约是 -40 到 30，这里我们先 grid on 让其自动适应，或者手动微调
ylim([-40, 30]); 
grid on; ax.GridAlpha = 0.15;

title('\fontname{Microsoft YaHei}IAP 各阶趋势分布图 (单轴)', 'FontSize', 16, 'FontWeight', 'bold');

% 图例
h_box  = patch(ax, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out  = scatter(ax, NaN, NaN, ScatterSize, UniColor, 'filled', 'MarkerFaceAlpha', 0.6);
legend(ax, [h_box, h_med, h_mean, h_out], {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 10, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');

% 保存
OutFile = fullfile(OutputDir, 'IAP_Trend_Boxplot_SingleAxis.jpg');
exportgraphics(gcf, OutFile, 'Resolution', 300);
fprintf('>> 图片已保存: %s\n', OutFile);

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
        scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
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
