function Plot_Cumulative_Sum(dataset_name, state, target_comp)
% 功能: 整合绘制1-8阶累加和的空间分布图和箱线图
% 特性: 科学防过曝色棒、高透轻量级星云散点、修复强制加载HSLA串台Bug、标题下方横排高级图例
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')
%   target_comp: 目标变量 ('TSLA' 或 'HSLA')，默认为 'TSLA'

clc; close all;
fprintf('========== 整合绘制1-8阶累加和的空间分布图和箱线图 ==========\n');

if nargin < 3
    target_comp = 'TSLA'; % 默认优先加载 TSLA
end

%% 1. 加载趋势和振幅结果文件
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', dataset_name, state));

if ~exist(trend_file, 'file'), error('❌ 找不到趋势结果文件: %s', trend_file); end
fprintf('>> 加载趋势数据: %s ...\n', trend_file);
load(trend_file);

amplitude_dir = 'D:\work\Task_Convergence\Amplitude_Results';
amplitude_file = fullfile(amplitude_dir, sprintf('%s_%s_Amplitude.mat', dataset_name, state));

if ~exist(amplitude_file, 'file'), error('❌ 找不到振幅结果文件: %s', amplitude_file); end
fprintf('>> 加载振幅数据: %s ...\n', amplitude_file);
load(amplitude_file);

%% 2. 准备数据 (🚨 核心修复：精准锁定变量，防止交叉感染)
fprintf('>> 🎯 当前锁定分析目标: 【%s】\n', target_comp);

if strcmpi(target_comp, 'TSLA') && exist('trend_TSLA', 'var') && exist('amp_TSLA', 'var')
    trend_data = trend_TSLA;
    sig_data = sig_TSLA;
    amplitude_data = amp_TSLA;
    comp_name = 'TSLA (Thermosteric)';
    fprintf('>> ✅ 成功提取纯温度 TSLA 的趋势、显著性与振幅数据。\n');
elseif strcmpi(target_comp, 'HSLA') && exist('trend_HSLA', 'var') && exist('amp_HSLA', 'var')
    trend_data = trend_HSLA;
    sig_data = sig_HSLA;
    amplitude_data = amp_HSLA;
    comp_name = 'HSLA (Halosteric)';
    fprintf('>> ✅ 成功提取盐收缩 HSLA 的趋势、显著性与振幅数据。\n');
else
    error('❌ 无法找到请求的变量 (%s)，或者趋势与振幅数据不完整！', target_comp);
end

% 经纬度处理
if ~exist('Lon', 'var') && exist('target_lon', 'var'), Lon = target_lon; end
if ~exist('Lat', 'var') && exist('target_lat', 'var'), Lat = target_lat; end
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end

% 确保Lon和Lat是向量
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% ⚠️ 第一步：先转置使得矩阵符合地图绘制要求 [Lat x Lon]
if size(trend_data, 1) ~= length(Lat)
    temp_t = nan(length(Lat), length(Lon), size(trend_data, 3));
    temp_s = nan(length(Lat), length(Lon), size(trend_data, 3));
    for i = 1:size(trend_data, 3)
        temp_t(:,:,i) = trend_data(:,:,i)';
        temp_s(:,:,i) = sig_data(:,:,i)';
    end
    trend_data = temp_t;
    sig_data = temp_s;
end

if size(amplitude_data, 1) ~= length(Lat)
    temp_a = nan(length(Lat), length(Lon), size(amplitude_data, 3));
    for i = 1:size(amplitude_data, 3)
        temp_a(:,:,i) = amplitude_data(:,:,i)';
    end
    amplitude_data = temp_a;
end

% ⚠️ 第二步：再进行经度 [0, 360] -> [-180, 180] 转换与切片拼接
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    for i = 1:size(trend_data, 3)
        trend_data(:,:,i) = trend_data(:, sort_idx, i);
        sig_data(:,:,i) = sig_data(:, sort_idx, i);
    end
    for i = 1:size(amplitude_data, 3)
        amplitude_data(:,:,i) = amplitude_data(:, sort_idx, i);
    end
end
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

%% 3. 计算1-8阶累加和
fprintf('>> 计算1-8阶累加和...\n');
max_terms = size(trend_data, 3);
if max_terms > 8, max_terms = 8; end

% 计算趋势累加和
trend_cumulative = zeros(size(trend_data, 1), size(trend_data, 2), max_terms);
trend_cumulative(:,:,1) = trend_data(:,:,1);
for i = 2:max_terms
    trend_cumulative(:,:,i) = trend_cumulative(:,:,i-1) + trend_data(:,:,i);
end

% 计算振幅累加和 (取绝对值后累加)
amplitude_cumulative = zeros(size(amplitude_data, 1), size(amplitude_data, 2), max_terms);
amplitude_cumulative(:,:,1) = abs(amplitude_data(:,:,1));
for i = 2:max_terms
    amplitude_cumulative(:,:,i) = amplitude_cumulative(:,:,i-1) + abs(amplitude_data(:,:,i));
end

%% 4. 绘制趋势累加和空间分布图 (顶刊 RdBu 色棒)
fprintf('>> 绘制趋势累加和空间分布图...\n');

figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Trend Cumulative Sum Spatial Maps');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = trend_cumulative(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % 显著性打点 (使用最后一阶的显著性)
    sig_mask = sig_data(:,:,max_terms); skip = 4; 
    [row, ~] = find(sig_mask(1:skip:end, 1:skip:end) == 1);
    if ~isempty(row)
        sig_lon = Lon_Grid(1:skip:end, 1:skip:end); sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==1), sig_lat(sig_mask(1:skip:end, 1:skip:end)==1), ...
               'o', 'markerfacecolor', [.3 .3 .3], 'markeredgecolor', 'none', 'markersize', 0.6);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap);
    
    % 🌟 核心优化 1：科学抗过曝截断
    if strcmp(state, 'Average'), pct = 99.0; else, pct = 99.5; end
    raw_max = prctile(abs(data(:)), pct);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    
    exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent; caxis([-final_limit, final_limit]); 
    
    title(sprintf('Cumulative Sum (1-%d Terms)', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) 1-%d Terms %s Trend Cumulative Sum', dataset_name, state, max_terms, strtok(comp_name, ' ')), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 5. 绘制振幅累加和空间分布图 (顶刊 振幅热力 色棒)
fprintf('>> 绘制振幅累加和空间分布图...\n');
figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Amplitude Cumulative Sum Spatial Maps');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊振幅专属色盘
Amp_Points = [1.00 1.00 1.00; 1.00 0.90 0.50; 0.95 0.50 0.15; 0.80 0.10 0.15; 0.40 0.00 0.00]; 
Amp_Locs = [0, 0.2, 0.5, 0.8, 1]; xx = linspace(0, 1, 256); 
MyAmpCmap = [interp1(Amp_Locs, Amp_Points(:,1), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,2), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = amplitude_cumulative(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyAmpCmap);
    
    if strcmp(state, 'Average'), pct = 99.0; else, pct = 99.5; end
    raw_max = prctile(data(:), pct); 
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    
    exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent; caxis([0, final_limit]); 
    
    title(sprintf('Cumulative Sum (1-%d Terms)', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm', 'FontSize', 8, 'FontWeight', 'bold'); cb.Ticks = linspace(0, final_limit, 5);
end

sgtitle(sprintf('%s (%s) 1-%d Terms %s Amplitude Cumulative Sum', dataset_name, state, max_terms, strtok(comp_name, ' ')), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 6. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(trend_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Trend_Cumulative = nan(num_valid, max_terms);
Data_Amplitude_Cumulative = nan(num_valid, max_terms);

for n = 1:max_terms
    data_slice_t = trend_cumulative(:,:,n); data_slice_a = amplitude_cumulative(:,:,n);
    vec_t = nan(num_valid, 1); vec_a = nan(num_valid, 1);
    for k = 1:num_valid
        vec_t(k) = data_slice_t(rows(k), cols(k));
        vec_a(k) = data_slice_a(rows(k), cols(k));
    end
    Data_Trend_Cumulative(:,n) = vec_t; Data_Amplitude_Cumulative(:,n) = vec_a;
end

if max(abs(Data_Trend_Cumulative(:))) > 500, Data_Trend_Cumulative = Data_Trend_Cumulative / 1000; end
if max(Data_Amplitude_Cumulative(:)) > 500, Data_Amplitude_Cumulative = Data_Amplitude_Cumulative / 1000; end

%% 7. 绘制趋势累加和箱线图
fprintf('>> 绘制趋势累加和箱线图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Trend Cumulative Sum Boxplot');

C1 = [0.902, 0.294, 0.208]; C2 = [0.302, 0.455, 0.690]; C3 = [0.000, 0.627, 0.408]; BoxBgColor = [0.96 0.96 0.96]; 
BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.20; ScatterSize  = 10;
X_Limit = [0.4, max_terms + 1.0]; 

% 🌟 调低主绘图区高度，为顶部的标题和图例预留空间
MainPos = [0.07, 0.12, 0.80, 0.73]; 

ax1 = axes('Position', MainPos); hold(ax1, 'on');
positions = 1:max_terms; colors = [C1; C2; C3; C1; C2; C3; C1; C2];

for k = 1:max_terms
    draw_jitter_unit(ax1, positions(k), Data_Trend_Cumulative(:,k), BoxWidth, ScatterWidth, ScatterSize, colors(mod(k-1,size(colors,1))+1,:), BoxBgColor, ScatterAlpha);
end

min_v = min(Data_Trend_Cumulative(:)); max_v = max(Data_Trend_Cumulative(:));
if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12, min_v = -1; max_v = 1; end
range = max_v - min_v; ylim(ax1, [min_v - range * 0.1, max_v + range * 0.1]); 
set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel(ax1, 'Trend Cumulative Sum (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);

xline(ax1, 0.7, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
set(ax1, 'XTick', positions, 'XTickLabel', strtrim(cellstr(num2str(positions')))', 'FontWeight', 'bold');
xlabel(ax1, 'Number of Terms', 'FontSize', 16, 'FontWeight', 'bold');

% ==========================================================
% 👑 标题与恢复后的横向排布顶级图例 (趋势箱线图)
% ==========================================================
ax_title = axes('Position', [0, 0.92, 1, 0.05], 'Visible', 'off');
text(ax_title, 0.5, 0.5, sprintf('%s (%s) %s Trend Cumulative Sum Distribution', dataset_name, state, strtok(comp_name, ' ')), 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');

ax_lgd = axes('Position', [MainPos(1), 0.87, MainPos(3), 0.04], 'XLim', [0, 100], 'YLim', [0, 1], 'Visible', 'off'); hold(ax_lgd, 'on');
patch(ax_lgd, [3, 5, 5, 3], [0.3, 0.3, 0.7, 0.7], BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
text(ax_lgd, 6, 0.5, 'IQR Box', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [19, 22], [0.5, 0.5], 'r-', 'LineWidth', 2);
text(ax_lgd, 23, 0.5, 'Median', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [34], [0.5], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
text(ax_lgd, 35, 0.5, 'Mean', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [45], [0.5], 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
text(ax_lgd, 46, 0.5, 'Max', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [55, 58], [0.5, 0.5], 'k-', 'LineWidth', 1.5);
text(ax_lgd, 59, 0.5, '2.5% \sim 97.5%', 'Interpreter', 'tex', 'FontSize', 13, 'FontName', 'Times New Roman');

scatter(ax_lgd, 79, 0.5, 50, C1, 'filled', 'MarkerFaceAlpha', 0.6);
text(ax_lgd, 81, 0.5, 'Outliers', 'FontSize', 13, 'FontName', 'Times New Roman');

hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on'; hPan = pan(gcf); hPan.Motion = 'both';

%% 8. 绘制振幅累加和箱线图
fprintf('>> 绘制振幅累加和箱线图...\n');
figure('Position', [50, 50, 1100, 500], 'Color', 'w', 'Name', 'Amplitude Cumulative Sum Boxplot');
ax2 = axes('Position', MainPos); hold(ax2, 'on');

for k = 1:max_terms
    draw_jitter_unit(ax2, positions(k), Data_Amplitude_Cumulative(:,k), BoxWidth, ScatterWidth, ScatterSize, colors(mod(k-1,size(colors,1))+1,:), BoxBgColor, ScatterAlpha);
end

min_v = min(Data_Amplitude_Cumulative(:)); max_v = max(Data_Amplitude_Cumulative(:));
if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12, min_v = 0; max_v = 1; end
range = max_v - min_v; ylim(ax2, [min_v - range * 0.1, max_v + range * 0.1]); 
set(ax2, 'XLim', X_Limit, 'YColor', C2, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel(ax2, 'Amplitude Cumulative Sum (mm)', 'FontWeight', 'bold', 'FontSize', 16);

xline(ax2, 0.7, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
set(ax2, 'XTick', positions, 'XTickLabel', strtrim(cellstr(num2str(positions')))', 'FontWeight', 'bold');
xlabel(ax2, 'Number of Terms', 'FontSize', 16, 'FontWeight', 'bold');

% ==========================================================
% 👑 标题与恢复后的横向排布顶级图例 (振幅箱线图)
% ==========================================================
ax_title2 = axes('Position', [0, 0.92, 1, 0.05], 'Visible', 'off');
text(ax_title2, 0.5, 0.5, sprintf('%s (%s) %s Amplitude Cumulative Sum Distribution', dataset_name, state, strtok(comp_name, ' ')), 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');

ax_lgd2 = axes('Position', [MainPos(1), 0.87, MainPos(3), 0.04], 'XLim', [0, 100], 'YLim', [0, 1], 'Visible', 'off'); hold(ax_lgd2, 'on');
patch(ax_lgd2, [3, 5, 5, 3], [0.3, 0.3, 0.7, 0.7], BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
text(ax_lgd2, 6, 0.5, 'IQR Box', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd2, [19, 22], [0.5, 0.5], 'r-', 'LineWidth', 2);
text(ax_lgd2, 23, 0.5, 'Median', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd2, [34], [0.5], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
text(ax_lgd2, 35, 0.5, 'Mean', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd2, [45], [0.5], 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
text(ax_lgd2, 46, 0.5, 'Max', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd2, [55, 58], [0.5, 0.5], 'k-', 'LineWidth', 1.5);
text(ax_lgd2, 59, 0.5, '2.5% \sim 97.5%', 'Interpreter', 'tex', 'FontSize', 13, 'FontName', 'Times New Roman');

% 振幅的散点通常以钴蓝 (C2) 为主基调展示
scatter(ax_lgd2, 79, 0.5, 50, C2, 'filled', 'MarkerFaceAlpha', 0.6);
text(ax_lgd2, 81, 0.5, 'Outliers', 'FontSize', 13, 'FontName', 'Times New Roman');

hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on'; hPan = pan(gcf); hPan.Motion = 'both';

%% 9. 保存结果图
fprintf('>> 保存结果图...\n');
save_dir = 'd:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i); if isvalid(fig)
        fig_name = get(fig, 'Name'); if isempty(fig_name), fig_name = sprintf('Figure_%d', i); end
        filename = sprintf('%s_%s_%s_%s.jpg', dataset_name, state, strtok(comp_name, ' '), strrep(fig_name, ' ', '_'));
        filepath = fullfile(save_dir, filename);
        saveas(fig, filepath, 'jpg'); fprintf('>> 已保存: %s\n', filepath);
    end
end
fprintf('>> 绘图完成！\n');
end

%% ============================================================
%% 内部函数
%% ============================================================
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    max_val = max(data); 
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    
    if ~isempty(outliers)
        if length(outliers) > 3000
            outliers = outliers(randperm(length(outliers), 3000));
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
    
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(ax, center, max_val, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6); 
end