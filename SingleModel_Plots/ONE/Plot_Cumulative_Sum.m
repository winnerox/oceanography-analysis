function Plot_Cumulative_Sum(dataset_name, state)
% 功能: 整合绘制1-8阶累加和的空间分布图和箱线图
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 整合绘制1-8阶累加和的空间分布图和箱线图 ==========\n');

%% 1. 加载趋势和振幅结果文件
% 加载趋势数据
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', dataset_name, state));

if ~exist(trend_file, 'file')
    error('❌ 找不到趋势结果文件: %s', trend_file);
end

fprintf('>> 加载趋势数据: %s ...\n', trend_file);
load(trend_file);

% 加载振幅数据
amplitude_dir = 'D:\work\Task_Convergence\Amplitude_Results';
amplitude_file = fullfile(amplitude_dir, sprintf('%s_%s_Amplitude.mat', dataset_name, state));

if ~exist(amplitude_file, 'file')
    error('❌ 找不到振幅结果文件: %s', amplitude_file);
end

fprintf('>> 加载振幅数据: %s ...\n', amplitude_file);
load(amplitude_file);

%% 2. 准备数据
% 提取趋势数据 (优先使用 HSLA 以显示真实变暖信号，若无则使用 TSLA)
if exist('trend_HSLA', 'var') && ~isempty(trend_HSLA)
    trend_data = trend_HSLA;
    sig_data = sig_HSLA;
    fprintf('>> 检测到 HSLA (总比容)，将绘制 HSLA 趋势以体现真实信号。\n');
else
    trend_data = trend_TSLA;
    sig_data = sig_TSLA;
    fprintf('>> 未检测到 HSLA，将绘制纯温度 TSLA 趋势。\n');
end

% 提取振幅数据 (优先使用 HSLA 以显示真实信号，若无则使用 TSLA)
if exist('amp_HSLA', 'var') && ~isempty(amp_HSLA)
    amplitude_data = amp_HSLA;
    fprintf('>> 检测到 HSLA (总比容)，将绘制 HSLA 振幅。\n');
else
    amplitude_data = amp_TSLA;
    fprintf('>> 未检测到 HSLA，将绘制纯温度 TSLA 振幅。\n');
end

% 经纬度处理
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end

% 确保Lon和Lat是向量
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180] 完美对齐
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

% 转置使得画图时不报错
if size(trend_data, 1) ~= length(Lat)
    for i = 1:size(trend_data, 3)
        temp_t(:,:,i) = trend_data(:,:,i)';
        temp_s(:,:,i) = sig_data(:,:,i)';
    end
    trend_data = temp_t;
    sig_data = temp_s;
end

if size(amplitude_data, 1) ~= length(Lat)
    for i = 1:size(amplitude_data, 3)
        temp_a(:,:,i) = amplitude_data(:,:,i)';
    end
    amplitude_data = temp_a;
end

%% 3. 计算1-8阶累加和
fprintf('>> 计算1-8阶累加和...\n');
max_terms = 8;

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
fprintf('>> 绘制全部八阶趋势累加和...\n');

% 顶刊发散型色盘 (深蓝-浅蓝-纯白-浅红-深红)
RdBu_ColorPoints = [
    0.1 0.2 0.5;
    0.6 0.8 0.9;
    1.0 1.0 1.0;
    1.0 0.7 0.7;
    0.6 0.1 0.1
];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; 
xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = trend_cumulative(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % 显著性打点 (使用最后一阶的显著性)
    sig_mask = sig_data(:,:,max_terms);
    skip = 3; 
    [row, col] = find(sig_mask(1:skip:end, 1:skip:end) == 1);
    if ~isempty(row)
        sig_lon = Lon_Grid(1:skip:end, 1:skip:end);
        sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==1), ...
               sig_lat(sig_mask(1:skip:end, 1:skip:end)==1), ...
               'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, MyDivergentCmap);
    
    raw_max = prctile(abs(data(:)), 99.8);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([-final_limit, final_limit]); 
    
    title(sprintf('Cumulative Sum (1-%d Terms)', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) 1-8 Terms Trend Cumulative Sum (p < 0.05)', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 5. 绘制振幅累加和空间分布图 (顶刊 振幅热力 色棒)
fprintf('>> 绘制振幅累加和空间分布图...\n');

figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Amplitude Cumulative Sum Spatial Maps');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊振幅专属色盘：白 -> 淡黄 -> 亮橘 -> 深红 
Amp_Points = [ 
    1.00, 1.00, 1.00;  % 0.0: 纯白 (背景/极弱区) 
    1.00, 0.90, 0.50;  % 0.2: 暖黄 (弱振幅) 
    0.95, 0.50, 0.15;  % 0.5: 亮橙 (中等振幅) 
    0.80, 0.10, 0.15;  % 0.8: 强红 (强振幅) 
    0.40, 0.00, 0.00   % 1.0: 暗红 (极值核心区) 
]; 
Amp_Locs = [0, 0.2, 0.5, 0.8, 1]; 
xx = linspace(0, 1, 256); 
MyAmpCmap = [interp1(Amp_Locs, Amp_Points(:,1), xx, 'pchip')', ... 
             interp1(Amp_Locs, Amp_Points(:,2), xx, 'pchip')', ... 
             interp1(Amp_Locs, Amp_Points(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = amplitude_cumulative(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, MyAmpCmap);
    
    raw_max = prctile(data(:), 99.5); % 振幅用 99.5% 截断更佳
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    
    caxis([0, final_limit]); 
    
    title(sprintf('Cumulative Sum (1-%d Terms)', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(0, final_limit, 5);
end

sgtitle(sprintf('%s (%s) 1-8 Terms Amplitude Cumulative Sum', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 6. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(trend_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Trend_Cumulative = nan(num_valid, 8);
Data_Amplitude_Cumulative = nan(num_valid, 8);

for n = 1:8
    % 趋势累加和数据
    data_slice_t = trend_cumulative(:,:,n);
    vec_t = nan(num_valid, 1);
    for k = 1:num_valid
        vec_t(k) = data_slice_t(rows(k), cols(k));
    end
    Data_Trend_Cumulative(:,n) = vec_t;
    
    % 振幅累加和数据
    data_slice_a = amplitude_cumulative(:,:,n);
    vec_a = nan(num_valid, 1);
    for k = 1:num_valid
        vec_a(k) = data_slice_a(rows(k), cols(k));
    end
    Data_Amplitude_Cumulative(:,n) = vec_a;
end

% 数据单位调整
if max(abs(Data_Trend_Cumulative(:))) > 500
     Data_Trend_Cumulative = Data_Trend_Cumulative / 1000;
end

if max(Data_Amplitude_Cumulative(:)) > 500
     Data_Amplitude_Cumulative = Data_Amplitude_Cumulative / 1000;
end

%% 7. 绘制趋势累加和箱线图
fprintf('>> 绘制趋势累加和箱线图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Trend Cumulative Sum Boxplot');

% 顶刊配色
C1 = [0.902, 0.294, 0.208]; % 樱桃红
C2 = [0.302, 0.455, 0.690]; % 钴蓝
C3 = [0.000, 0.627, 0.408]; % 薄荷绿
BoxBgColor = [0.96 0.96 0.96]; % 极浅高级灰底色

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

% 设置箱线图位置和范围
X_Limit = [0.4, 9.0];
MainPos = [0.07, 0.14, 0.80, 0.80]; 

% 绘制箱线图 (使用单个y轴，因为累加和的量级差异不大)
ax1 = axes('Position', MainPos); hold(ax1, 'on');

% 为每个累加和绘制箱线图
positions = 1:8;
colors = [C1; C2; C3; C1; C2; C3; C1; C2];

for k = 1:8
    draw_jitter_unit(ax1, positions(k), Data_Trend_Cumulative(:,k), BoxWidth, ScatterWidth, ScatterSize, colors(k,:), BoxBgColor, ScatterAlpha);
end

% 使用数据的实际最大最小值
min_v = min(Data_Trend_Cumulative(:));
max_v = max(Data_Trend_Cumulative(:));
if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12
    min_v = -1;
    max_v = 1;
end
% 添加10%的边距
range = max_v - min_v;
Limit_min = min_v - range * 0.1;
Limit_max = max_v + range * 0.1;
ylim(ax1, [Limit_min, Limit_max]); 
set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel(ax1, 'Trend Cumulative Sum (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);

% --- 装饰与分割线 ---
xline(ax1, 0.7, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

All_Ticks = positions;
All_Labels = {'1','2','3', '4','5', '6','7','8'};
set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
xlabel(ax1, 'Number of Terms', 'FontSize', 16, 'FontWeight', 'bold');
title(ax1, sprintf('%s (%s) Trend Cumulative Sum Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');

% === 图例 ===
h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out  = scatter(ax1, NaN, NaN, ScatterSize, C1, 'filled', 'MarkerFaceAlpha', 0.6);

legend(ax1, [h_box, h_med, h_mean, h_out], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'none');

% === 交互优化 ===
hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
hPan = pan(gcf); hPan.Motion = 'both';

%% 8. 绘制振幅累加和箱线图
fprintf('>> 绘制振幅累加和箱线图...\n');
figure('Position', [50, 50, 1100, 500], 'Color', 'w', 'Name', 'Amplitude Cumulative Sum Boxplot');

% 绘制箱线图 (使用单个y轴)
ax2 = axes('Position', MainPos); hold(ax2, 'on');

for k = 1:8
    draw_jitter_unit(ax2, positions(k), Data_Amplitude_Cumulative(:,k), BoxWidth, ScatterWidth, ScatterSize, colors(k,:), BoxBgColor, ScatterAlpha);
end

% 使用数据的实际最大最小值
min_v = min(Data_Amplitude_Cumulative(:));
max_v = max(Data_Amplitude_Cumulative(:));
if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12
    min_v = 0;
    max_v = 1;
end
% 添加10%的边距
range = max_v - min_v;
Limit_min = min_v - range * 0.1;
Limit_max = max_v + range * 0.1;
ylim(ax2, [Limit_min, Limit_max]); 
set(ax2, 'XLim', X_Limit, 'YColor', C2, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel(ax2, 'Amplitude Cumulative Sum (mm)', 'FontWeight', 'bold', 'FontSize', 16);

% --- 装饰与分割线 ---
xline(ax2, 0.7, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

set(ax2, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
xlabel(ax2, 'Number of Terms', 'FontSize', 16, 'FontWeight', 'bold');
title(ax2, sprintf('%s (%s) Amplitude Cumulative Sum Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');

% === 图例 ===
h_box  = patch(ax2, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax2, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax2, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out  = scatter(ax2, NaN, NaN, ScatterSize, C2, 'filled', 'MarkerFaceAlpha', 0.6);

legend(ax2, [h_box, h_med, h_mean, h_out], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 11, 'Box', 'on', 'EdgeColor', 'none');

% === 交互优化 ===
hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
hPan = pan(gcf); hPan.Motion = 'both';

%% 9. 保存结果图
fprintf('>> 保存结果图...\n');

% 创建保存目录
save_dir = 'd:\work\Figures';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('>> 创建保存目录: %s\n', save_dir);
end

% 保存所有打开的图形窗口
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i);
    if isvalid(fig)
        % 获取图形名称
        fig_name = get(fig, 'Name');
        if isempty(fig_name)
            fig_name = sprintf('Figure_%d', i);
        end
        
        % 构建文件名
        filename = sprintf('%s_%s_%s.jpg', dataset_name, state, strrep(fig_name, ' ', '_'));
        filepath = fullfile(save_dir, filename);
        
        % 保存图片
        saveas(fig, filepath, 'jpg');
        fprintf('>> 已保存: %s\n', filepath);
    end
end

fprintf('>> 绘图完成！\n');
fprintf('>> 图片已保存到: %s\n', save_dir);
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
    max_val = max(data); % 计算最大值
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    
    % 绘制散点 (去除黑边框，恢复果冻质感透明叠加)
    if ~isempty(outliers)
        if length(outliers) > 2000
            outliers = outliers(randperm(length(outliers), 2000));
        end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        try
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s);
        end
    end
    
    % 绘制须线 (恢复纤细的 1.2 细线)
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    % 绘制箱体
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, ...
        'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    % 绘制中位线、均值点和最大值点
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(ax, center, max_val, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % 绘制最大值点
end
