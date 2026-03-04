function Plot_Trend_Combined(dataset_name, state)
% 功能: 整合绘制趋势空间分布图和箱线图
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 整合绘制趋势空间分布图和箱线图 ==========\n');

%% 1. 加载趋势结果文件
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', dataset_name, state));

if ~exist(trend_file, 'file')
    error('❌ 找不到趋势结果文件: %s', trend_file);
end

fprintf('>> 加载趋势数据: %s ...\n', trend_file);
load(trend_file);

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

%% 3. 绘制空间分布图 (顶刊 RdBu 色棒)
fprintf('>> 绘制趋势空间分布图...\n');

% 确定显示的阶数
max_terms = 8;
figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Spatial Maps');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
fprintf('>> 绘制全部八阶趋势...\n');

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
    
    data = trend_data(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % 显著性打点
    sig_mask = sig_data(:,:,term_idx);
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
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) 1-8 Terms Trend Distribution (p < 0.05)', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 4. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(trend_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Trend = nan(num_valid, 8);

for n = 1:8
    data_slice = trend_data(:,:,n);
    vec_t = nan(num_valid, 1);
    for k = 1:num_valid
        vec_t(k) = data_slice(rows(k), cols(k));
    end
    Data_Trend(:,n) = vec_t;
end

if max(abs(Data_Trend(:))) > 500
     Data_Trend = Data_Trend / 1000;
end

%% 5. 绘制箱线图
fprintf('>> 绘制趋势箱线图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Trend Boxplot');

% 顶刊配色
C1 = [0.902, 0.294, 0.208]; % 樱桃红
C2 = [0.302, 0.455, 0.690]; % 钴蓝
C3 = [0.000, 0.627, 0.408]; % 薄荷绿
BoxBgColor = [0.96 0.96 0.96]; % 极浅高级灰底色

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

% 根据状态设置不同的划分逻辑
X_Limit = [0.4, 9.0];
MainPos = [0.07, 0.14, 0.80, 0.80]; 

if strcmp(state, 'Average')
    % 平均态划分：按量级分组
    Pos_G1 = [1];      % T1 (主趋势，量级 ~10^1)
    Pos_G2 = [2, 3];   % T2-T3 (量级 ~10^0 到 10^-1)
    Pos_G3 = [4.2, 5.2]; % T4-T5 (量级 ~10^-3)
    Pos_G4 = [6.4, 7.4, 8.4]; % T6-T8 (量级 ~10^-5 到 10^-7)
else
    % 标准态划分：按阶数分组
    Pos_G1 = [1, 2, 3];   % T1-T3 (前三个主项)
    Pos_G2 = [4.2, 5.2]; % T4-T5 (中间两项)
    Pos_G3 = [6.4, 7.4, 8.4]; % T6-T8 (最后三项)
end 

if strcmp(state, 'Average')
    % 平均态绘制逻辑：四个y轴
    % --- Layer 1 (Red) - T1 主趋势 ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    draw_jitter_unit(ax1, Pos_G1(1), Data_Trend(:,1), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
    % 使用数据的实际最大最小值
    min_v1 = min(Data_Trend(:,1));
    max_v1 = max(Data_Trend(:,1));
    if isnan(min_v1)||isnan(max_v1)||(max_v1-min_v1)<1e-12
        min_v1 = -1;
        max_v1 = 1;
    end
    % 添加10%的边距
    range1 = max_v1 - min_v1;
    Limit1_min = min_v1 - range1 * 0.1;
    Limit1_max = max_v1 + range1 * 0.1;
    ylim(ax1, [Limit1_min, Limit1_max]); 
    set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax1, 'Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 2 (Blue) - T2-T3 趋势 ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        draw_jitter_unit(ax2, Pos_G2(k), Data_Trend(:,k+1), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v2 = min(Data_Trend(:,2:3), [], 'all');
    max_v2 = max(Data_Trend(:,2:3), [], 'all');
    if isnan(min_v2)||isnan(max_v2)||(max_v2-min_v2)<1e-12
        min_v2 = -1;
        max_v2 = 1;
    end
    % 添加10%的边距
    range2 = max_v2 - min_v2;
    Limit2_min = min_v2 - range2 * 0.1;
    Limit2_max = max_v2 + range2 * 0.1;
    ylim(ax2, [Limit2_min, Limit2_max]); 
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- Layer 3 (Green) - T4-T5 趋势 ---
    Offset_Ratio = 0.04; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = X_Limit * Scale_Factor; 
    
    for k = 1:2
        draw_jitter_unit(ax3, Pos_G3(k), Data_Trend(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v3 = min(Data_Trend(:,4:5), [], 'all');
    max_v3 = max(Data_Trend(:,4:5), [], 'all');
    if isnan(min_v3)||isnan(max_v3)||(max_v3-min_v3)<1e-12
        min_v3 = -1;
        max_v3 = 1;
    end
    % 添加10%的边距
    range3 = max_v3 - min_v3;
    Limit3_min = min_v3 - range3 * 0.1;
    Limit3_max = max_v3 + range3 * 0.1;
    ylim(ax3, [Limit3_min, Limit3_max]);
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- Layer 4 (Purple) - T6-T8 趋势 ---
    Offset_Ratio2 = 0.08; 
    Pos4 = MainPos; Pos4(3) = MainPos(3) * (1 + Offset_Ratio2);
    ax4 = axes('Position', Pos4); hold(ax4, 'on');
    Scale_Factor2 = Pos4(3) / MainPos(3); 
    Real_X_Limit2 = X_Limit * Scale_Factor2; 
    
    for k = 1:3
        draw_jitter_unit(ax4, Pos_G4(k), Data_Trend(:,k+5), BoxWidth, ScatterWidth, ScatterSize, [0.5 0 0.5], BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v4 = min(Data_Trend(:,6:8), [], 'all');
    max_v4 = max(Data_Trend(:,6:8), [], 'all');
    if isnan(min_v4)||isnan(max_v4)||(max_v4-min_v4)<1e-12
        min_v4 = -1;
        max_v4 = 1;
    end
    % 添加10%的边距
    range4 = max_v4 - min_v4;
    Limit4_min = min_v4 - range4 * 0.1;
    Limit4_max = max_v4 + range4 * 0.1;
    ylim(ax4, [Limit4_min, Limit4_max]);
    set(ax4, 'Color', 'none', 'XLim', Real_X_Limit2, 'YAxisLocation', 'right', 'YColor', [0.5 0 0.5], 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- 装饰与分割线 ---
    xline(ax1, 1.5, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    
    All_Ticks = [Pos_G1, Pos_G2, Pos_G3, Pos_G4];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    title(ax1, sprintf('%s (%s) Terms Trend Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');
    
    % === 图例 ===
    h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
    h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
    h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    h_out  = scatter(ax1, NaN, NaN, ScatterSize, C1, 'filled', 'MarkerFaceAlpha', 0.6);
    
    legend(ax1, [h_box, h_med, h_mean, h_out], ...
           {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
           'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'none');
    
    % === 交互优化 ===
    % 四个y轴，需要同步
    % 计算比例因子
    range1 = Limit1_max - Limit1_min;
    range2 = Limit2_max - Limit2_min;
    range3 = Limit3_max - Limit3_min;
    range4 = Limit4_max - Limit4_min;
    Ratio2 = range2 / range1;
    Ratio3 = range3 / range1;
    Ratio4 = range4 / range1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncAllAxes4(ax1, ax2, ax3, ax4, Ratio2, Ratio3, Ratio4, 1, Scale_Factor, Scale_Factor2);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncAllAxes4(ax1, ax2, ax3, ax4, Ratio2, Ratio3, Ratio4, 1, Scale_Factor, Scale_Factor2);
else
    % 标准态绘制逻辑：三个y轴
    % --- Layer 1 (Red) - T1-T3 趋势 ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for k = 1:3
        draw_jitter_unit(ax1, Pos_G1(k), Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v1 = min(Data_Trend(:,1:3), [], 'all');
    max_v1 = max(Data_Trend(:,1:3), [], 'all');
    if isnan(min_v1)||isnan(max_v1)||(max_v1-min_v1)<1e-12
        min_v1 = -1;
        max_v1 = 1;
    end
    % 添加10%的边距
    range1 = max_v1 - min_v1;
    Limit1_min = min_v1 - range1 * 0.1;
    Limit1_max = max_v1 + range1 * 0.1;
    ylim(ax1, [Limit1_min, Limit1_max]); 
    set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax1, 'Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 2 (Blue) - T4-T5 趋势 ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        draw_jitter_unit(ax2, Pos_G2(k), Data_Trend(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v2 = min(Data_Trend(:,4:5), [], 'all');
    max_v2 = max(Data_Trend(:,4:5), [], 'all');
    if isnan(min_v2)||isnan(max_v2)||(max_v2-min_v2)<1e-12
        min_v2 = -1;
        max_v2 = 1;
    end
    % 添加10%的边距
    range2 = max_v2 - min_v2;
    Limit2_min = min_v2 - range2 * 0.1;
    Limit2_max = max_v2 + range2 * 0.1;
    ylim(ax2, [Limit2_min, Limit2_max]); 
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- Layer 3 (Green) - T6-T8 趋势 ---
    Offset_Ratio = 0.04; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = X_Limit * Scale_Factor; 
    
    for k = 1:3
        draw_jitter_unit(ax3, Pos_G3(k), Data_Trend(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
    end
    % 使用数据的实际最大最小值
    min_v3 = min(Data_Trend(:,6:8), [], 'all');
    max_v3 = max(Data_Trend(:,6:8), [], 'all');
    if isnan(min_v3)||isnan(max_v3)||(max_v3-min_v3)<1e-12
        min_v3 = -1;
        max_v3 = 1;
    end
    % 添加10%的边距
    range3 = max_v3 - min_v3;
    Limit3_min = min_v3 - range3 * 0.1;
    Limit3_max = max_v3 + range3 * 0.1;
    ylim(ax3, [Limit3_min, Limit3_max]);
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- 装饰与分割线 ---
    xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    
    All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    title(ax1, sprintf('%s (%s) Terms Trend Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');
    
    % === 图例 ===
    h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
    h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
    h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    h_out  = scatter(ax1, NaN, NaN, ScatterSize, C1, 'filled', 'MarkerFaceAlpha', 0.6);
    
    legend(ax1, [h_box, h_med, h_mean, h_out], ...
           {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
           'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'none');
    
    % === 交互优化 ===
    % 三个y轴，需要同步
    % 计算比例因子
    range1 = Limit1_max - Limit1_min;
    range2 = Limit2_max - Limit2_min;
    range3 = Limit3_max - Limit3_min;
    Ratio2 = range2 / range1;
    Ratio3 = range3 / range1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncThreeAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncThreeAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
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
    
    % 绘制中位线和均值点
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

function SyncAllAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        % 获取当前 ax1 的 Y 轴范围
        y1_lim = get(ax1, 'YLim');
        y1_min = y1_lim(1);
        y1_max = y1_lim(2);
        
        % 计算 ax2 和 ax3 的 Y 轴范围
        % 保持比例关系
        y2_min = y1_min * R2;
        y2_max = y1_max * R2;
        y3_min = y1_min * R3;
        y3_max = y1_max * R3;
        
        % 设置 Y 轴范围
        set(ax2, 'YLim', [y2_min, y2_max]);
        set(ax3, 'YLim', [y3_min, y3_max]);
        set(ax3, 'XLim', get(ax1, 'XLim') * Factor);
    catch ME
        % 捕获并显示错误，以便调试
        disp(['SyncAllAxes 错误: ', ME.message]);
    end
end

function SyncAllAxes4(ax1, ax2, ax3, ax4, R2, R3, R4, Factor1, Factor2)
    try
        % 获取当前 ax1 的 Y 轴范围
        y1_lim = get(ax1, 'YLim');
        y1_min = y1_lim(1);
        y1_max = y1_lim(2);
        
        % 计算 ax2、ax3 和 ax4 的 Y 轴范围
        % 保持比例关系
        y2_min = y1_min * R2;
        y2_max = y1_max * R2;
        y3_min = y1_min * R3;
        y3_max = y1_max * R3;
        y4_min = y1_min * R4;
        y4_max = y1_max * R4;
        
        % 设置 Y 轴范围
        set(ax2, 'YLim', [y2_min, y2_max]);
        set(ax3, 'YLim', [y3_min, y3_max]);
        set(ax4, 'YLim', [y4_min, y4_max]);
        set(ax3, 'XLim', get(ax1, 'XLim') * Factor1);
        set(ax4, 'XLim', get(ax1, 'XLim') * Factor2);
    catch ME
        % 捕获并显示错误，以便调试
        disp(['SyncAllAxes4 错误: ', ME.message]);
    end
end

function SyncTwoAxes(ax1, ax2, R2, Factor)
    try
        % 获取当前 ax1 的 Y 轴范围
        y1_lim = get(ax1, 'YLim');
        y1_min = y1_lim(1);
        y1_max = y1_lim(2);
        
        % 计算 ax2 的 Y 轴范围
        y2_min = y1_min * R2;
        y2_max = y1_max * R2;
        
        % 设置 Y 轴范围
        set(ax2, 'YLim', [y2_min, y2_max]);
        set(ax2, 'XLim', get(ax1, 'XLim') * Factor);
    catch ME
        disp(['SyncTwoAxes 错误: ', ME.message]);
    end
end

function SyncThreeAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        % 获取当前 ax1 的 Y 轴范围
        y1_lim = get(ax1, 'YLim');
        y1_min = y1_lim(1);
        y1_max = y1_lim(2);
        
        % 计算 ax2 和 ax3 的 Y 轴范围
        y2_min = y1_min * R2;
        y2_max = y1_max * R2;
        y3_min = y1_min * R3;
        y3_max = y1_max * R3;
        
        % 设置 Y 轴范围
        set(ax2, 'YLim', [y2_min, y2_max]);
        set(ax3, 'YLim', [y3_min, y3_max]);
        set(ax3, 'XLim', get(ax1, 'XLim') * Factor);
    catch ME
        disp(['SyncThreeAxes 错误: ', ME.message]);
    end
end

function SyncFourAxes(ax1, ax2, ax3, ax4, R2, R3, R4, Factor1, Factor2, Factor3)
    try
        % 获取当前 ax1 的 Y 轴范围
        y1_lim = get(ax1, 'YLim');
        y1_min = y1_lim(1);
        y1_max = y1_lim(2);
        
        % 计算 ax2、ax3 和 ax4 的 Y 轴范围
        y2_min = y1_min * R2;
        y2_max = y1_max * R2;
        y3_min = y1_min * R3;
        y3_max = y1_max * R3;
        y4_min = y1_min * R4;
        y4_max = y1_max * R4;
        
        % 设置 Y 轴范围
        set(ax2, 'YLim', [y2_min, y2_max]);
        set(ax3, 'YLim', [y3_min, y3_max]);
        set(ax4, 'YLim', [y4_min, y4_max]);
        set(ax2, 'XLim', get(ax1, 'XLim') * Factor1);
        set(ax3, 'XLim', get(ax1, 'XLim') * Factor2);
        set(ax4, 'XLim', get(ax1, 'XLim') * Factor3);
    catch ME
        disp(['SyncFourAxes 错误: ', ME.message]);
    end
end