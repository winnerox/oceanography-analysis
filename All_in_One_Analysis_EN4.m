%% All_in_One_Analysis_EN4_Fixed_v2.m
% =========================================================================
% 全能整合版 (v2 修复维度错误)：TSLA 1-8 阶 统计与绘图
% 修复日志：
%   1. [Fix] 修正 temp_trend 初始化维度为 [nlat, nlon]，解决赋值报错。
%   2. [Check] 再次确认单位为 mm，不重复乘 1000。
% =========================================================================

clear; clc; close all;

%% ========================================================================
%% 1. 数据加载与预处理
%% ========================================================================
InputFile = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms1to8_Average.mat';

if ~exist(InputFile, 'file')
    error('❌ 找不到文件: %s \n请检查路径是否正确。', InputFile);
end
fprintf('>> [1/4] 加载数据: %s ...\n', InputFile);
load(InputFile);

% 经纬度修正
if max(lon) > 180
    lon(lon > 180) = lon(lon > 180) - 360;
end
[lon, sort_idx] = sort(lon);
TSLA_AllOrders = TSLA_AllOrders(sort_idx, :, :, :); 

[nlon, nlat, ntime, max_order] = size(TSLA_AllOrders);
% t_vec 需要是列向量
t_vec = time_axis(:);

%% ========================================================================
%% 2. 核心计算 (Trend, Sig, Amp)
%% ========================================================================
fprintf('>> [2/4]正在并行计算 趋势(Trend) 和 振幅(Amp)...\n');

% 最终结果矩阵：[Lat, Lon, Order] -> [173, 360, 8]
Map_Trend = nan(nlat, nlon, max_order);
Map_Sig   = nan(nlat, nlon, max_order);
Map_Amp   = nan(nlat, nlon, max_order);

Box_Trend = cell(1, max_order);
Box_Amp   = cell(1, max_order);

for n = 1:max_order
    fprintf('   正在处理 Term %d / %d ...\n', n, max_order);
    
    slice = squeeze(TSLA_AllOrders(:,:,:,n)); % [Lon, Lat, Time]
    
    % 【关键修复】初始化维度必须是 [nlat, nlon] 即 [173, 360]
    temp_trend = nan(nlat, nlon);
    temp_sig   = nan(nlat, nlon);
    temp_amp   = nan(nlat, nlon);
    
    for i = 1:nlon
        for j = 1:nlat
            ts = squeeze(slice(i,j,:));
            idx_valid = ~isnan(ts);
            
            if sum(idx_valid) < length(ts) * 0.7, continue; end
            
            % --- 计算趋势 ---
            [p_coeff, ~] = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
            trend_val = p_coeff(1); % 单位已是 mm/yr
            
            % --- 计算振幅 ---
            y_fit = polyval(p_coeff, t_vec(idx_valid));
            resid = ts(idx_valid) - y_fit;
            amp_val = std(resid);   % 单位已是 mm
            
            % --- 显著性 ---
            n_val = sum(idx_valid);
            X = [ones(n_val,1) t_vec(idx_valid)];
            XTX_inv = inv(X' * X);
            se_slope = sqrt(var(resid) * XTX_inv(2,2));
            t_stat = p_coeff(1) / se_slope;
            p_val = 2 * (1 - tcdf(abs(t_stat), n_val - 2));
            
            % 存储：注意这里进行了转置 (j, i) -> (Lat, Lon)
            temp_trend(j,i) = trend_val; 
            temp_amp(j,i)   = amp_val;   
            if p_val < 0.05
                temp_sig(j,i) = 1;       
            end
        end
    end
    
    % 赋值回主矩阵 (维度匹配：173x360 = 173x360)
    Map_Trend(:,:,n) = temp_trend;
    Map_Sig(:,:,n)   = temp_sig;
    Map_Amp(:,:,n)   = temp_amp;
    
    vec_t = temp_trend(:); Box_Trend{n} = vec_t(~isnan(vec_t));
    vec_a = temp_amp(:);   Box_Amp{n}   = vec_a(~isnan(vec_a));
end

% 准备绘图网格 [Lat, Lon]
[Lon_Plot, Lat_Plot] = meshgrid(lon, lat); 

%% ========================================================================
%% 3. 绘图配置参数
%% ========================================================================
C1 = [0.85 0.325 0.098]; C2 = [0 0.447 0.741]; C3 = [0.466 0.674 0.188];
UniColor = [0.2 0.4 0.7]; BoxBgColor = [0.95 0.95 1];
BoxWidth = 0.55; ScatterWidth = 0.4; ScatterSize = 12; ScatterAlpha = 0.4;
Pos_G1 = [1, 2, 3]; Pos_G2 = [4.2, 5.2]; Pos_G3 = [6.4, 7.4, 8.4];
X_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};

%% ========================================================================
%% 4. 绘图 Part A: 振幅箱线图 (Amplitude)
%% ========================================================================
fprintf('>> [3/4] 绘制箱线图...\n');

% Figure 1: Amp Compact
figure('Position', [50, 500, 1200, 550], 'Color', 'w', 'Name', 'Amp Compact');
MainPos = [0.07, 0.14, 0.80, 0.80]; 
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:3, draw_jitter_unit(ax1, Pos_G1(k), Box_Amp{k}, BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha); end
max_v1 = prctile(cell2mat(Box_Amp(1:3)'), 99.8); ylim(ax1, [0, max_v1*1.1]);
set(ax1, 'XLim', [0.4 9.0], 'YColor', C1, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 1:2, draw_jitter_unit(ax2, Pos_G2(k), Box_Amp{k+3}, BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha); end
max_v2 = prctile(cell2mat(Box_Amp(4:5)'), 99.8); ylim(ax2, [0, max_v2*1.2]);
set(ax2, 'Color', 'none', 'XLim', [0.4 9.0], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

Pos3 = MainPos; Pos3(3) = MainPos(3) * 1.08; 
ax3 = axes('Position', Pos3); hold(ax3, 'on');
for k = 1:3, draw_jitter_unit(ax3, Pos_G3(k), Box_Amp{k+5}, BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha); end
max_v3 = prctile(cell2mat(Box_Amp(6:8)'), 99.8); ylim(ax3, [0, max_v3*1.2]);
set(ax3, 'Color', 'none', 'XLim', [0.4 9.0]*1.08, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

yline(ax1, 0, '-', 'Color', 'k'); 
xline(ax1, 3.6, ':', 'Color', [.5 .5 .5]); xline(ax1, 5.8, ':', 'Color', [.5 .5 .5]);
set(ax1, 'XTick', [Pos_G1, Pos_G2, Pos_G3], 'XTickLabel', X_Labels);
title(ax1, '\fontname{Microsoft YaHei}各阶独立项振幅 (紧凑布局)', 'FontSize', 16, 'FontWeight', 'bold');
add_legend(ax1, BoxBgColor, ScatterSize, 'Outliers (Colored)');

% Figure 2: Amp Single
figure('Position', [100, 100, 1200, 550], 'Color', 'w', 'Name', 'Amp Single');
ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');
for i = 1:8, draw_jitter_unit(ax, i, Box_Amp{i}, 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha); end
yline(ax, 0, '-', 'Color', 'k'); set(ax, 'XTick', 1:8, 'XTickLabel', X_Labels, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax, 'Amplitude (mm)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax, '\fontname{Microsoft YaHei}各阶独立项振幅 (单轴)', 'FontSize', 16, 'FontWeight', 'bold');
grid on; ax.GridAlpha = 0.15; add_legend(ax, BoxBgColor, ScatterSize, 'Outliers');

%% ========================================================================
%% 5. 绘图 Part B: 趋势箱线图 (Trend)
%% ========================================================================

% Figure 3: Trend Compact
figure('Position', [150, 500, 1200, 550], 'Color', 'w', 'Name', 'Trend Compact');
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:3, draw_jitter_unit(ax1, Pos_G1(k), Box_Trend{k}, BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha); end
max_v1 = prctile(abs(cell2mat(Box_Trend(1:3)')), 99.8); ylim(ax1, [-max_v1*1.1, max_v1*1.1]);
set(ax1, 'XLim', [0.4 9.0], 'YColor', C1, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 1:2, draw_jitter_unit(ax2, Pos_G2(k), Box_Trend{k+3}, BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha); end
max_v2 = prctile(abs(cell2mat(Box_Trend(4:5)')), 99.8); ylim(ax2, [-max_v2*1.2, max_v2*1.2]);
set(ax2, 'Color', 'none', 'XLim', [0.4 9.0], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

ax3 = axes('Position', Pos3); hold(ax3, 'on');
for k = 1:3, draw_jitter_unit(ax3, Pos_G3(k), Box_Trend{k+5}, BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha); end
max_v3 = prctile(abs(cell2mat(Box_Trend(6:8)')), 99.8); ylim(ax3, [-max_v3*1.2, max_v3*1.2]);
set(ax3, 'Color', 'none', 'XLim', [0.4 9.0]*1.08, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

yline(ax1, 0, '--', 'Color', [.6 .6 .6]); 
xline(ax1, 3.6, ':', 'Color', [.5 .5 .5]); xline(ax1, 5.8, ':', 'Color', [.5 .5 .5]);
set(ax1, 'XTick', [Pos_G1, Pos_G2, Pos_G3], 'XTickLabel', X_Labels);
title(ax1, '\fontname{Microsoft YaHei}各阶独立项趋势 (紧凑布局)', 'FontSize', 16, 'FontWeight', 'bold');
add_legend(ax1, BoxBgColor, ScatterSize, 'Outliers (Colored)');

% Figure 4: Trend Single
figure('Position', [200, 100, 1200, 550], 'Color', 'w', 'Name', 'Trend Single');
ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');
for i = 1:8, draw_jitter_unit(ax, i, Box_Trend{i}, 0.6, ScatterWidth, ScatterSize, UniColor, BoxBgColor, ScatterAlpha); end
yline(ax, 0, '--', 'Color', [.5 .5 .5]); set(ax, 'XTick', 1:8, 'XTickLabel', X_Labels, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax, 'Trend (mm/yr)', 'FontSize', 14, 'FontWeight', 'bold');
title(ax, '\fontname{Microsoft YaHei}各阶独立项趋势 (单轴)', 'FontSize', 16, 'FontWeight', 'bold');
grid on; ax.GridAlpha = 0.15; add_legend(ax, BoxBgColor, ScatterSize, 'Outliers');

%% ========================================================================
%% 6. 绘图 Part C: 空间分布图 (Spatial Maps)
%% ========================================================================
fprintf('>> [4/4] 绘制空间分布图...\n');

% 颜色表
MyBlue=[0.05 0.25 0.6]; MyLBlue=[0.6 0.8 1]; MyWhite=[1 1 1]; MyLRed=[1 0.7 0.6]; MyRed=[0.7 0.05 0.05];
TrendCmap = [interp1([0 .25 .5 .75 1], [MyBlue(:,1);MyLBlue(:,1);MyWhite(:,1);MyLRed(:,1);MyRed(:,1)], linspace(0,1,256),'pchip')', ...
             interp1([0 .25 .5 .75 1], [MyBlue(:,2);MyLBlue(:,2);MyWhite(:,2);MyLRed(:,2);MyRed(:,2)], linspace(0,1,256),'pchip')', ...
             interp1([0 .25 .5 .75 1], [MyBlue(:,3);MyLBlue(:,3);MyWhite(:,3);MyLRed(:,3);MyRed(:,3)], linspace(0,1,256),'pchip')'];
AmpRed = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
AmpCmap = [interp1([0 0.1 0.4 0.8 1], AmpRed(:,1), linspace(0,1,256))', interp1([0 0.1 0.4 0.8 1], AmpRed(:,2), linspace(0,1,256))', interp1([0 0.1 0.4 0.8 1], AmpRed(:,3), linspace(0,1,256))'];

% Figure 5: Map Amp
figure('Position', [250, 50, 1200, 1400], 'Color', 'w', 'Name', 'Map Amp');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for term_idx = 1:8
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Map_Amp(:,:,term_idx);
    m_pcolor(Lon_Plot, Lat_Plot, data); shading flat; 
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, AmpCmap);
    raw_max = prctile(data(:), 95.0); if raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); nice_base = ceil(raw_max/10^exponent*2)/2; 
    limit = nice_base * 10^exponent;
    caxis([0, limit]);
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; cb.FontSize=8; title(cb, 'mm', 'FontSize', 8); cb.Ticks = linspace(0, limit, 5);
end
sgtitle('\fontname{Microsoft YaHei}1-8 阶独立项振幅分布', 'FontSize', 16, 'FontWeight', 'bold');

% Figure 6: Map Trend
figure('Position', [300, 50, 1200, 1400], 'Color', 'w', 'Name', 'Map Trend');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for term_idx = 1:8
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Map_Trend(:,:,term_idx);
    m_pcolor(Lon_Plot, Lat_Plot, data); shading flat; hold on;
    
    % 打点
    sig_mask = Map_Sig(:,:,term_idx); skip=3;
    [row, col] = find(sig_mask(1:skip:end, 1:skip:end) == 1);
    if ~isempty(row)
         sig_lon = Lon_Plot(1:skip:end, 1:skip:end); sig_lat = Lat_Plot(1:skip:end, 1:skip:end);
         m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==1), sig_lat(sig_mask(1:skip:end, 1:skip:end)==1), ...
             'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
    
    colormap(gca, TrendCmap);
    raw_max = prctile(abs(data(:)), 99.8); if raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); nice_base = ceil(raw_max/10^exponent);
    limit = nice_base * 10^exponent;
    caxis([-limit, limit]);
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; cb.FontSize=8; title(cb, 'mm/yr', 'FontSize', 8); cb.Ticks = linspace(-limit, limit, 5);
end
sgtitle('\fontname{Microsoft YaHei}1-8 阶独立项趋势分布 (p < 0.05)', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> ✅ 所有绘图完成！(维度修复版)\n');

%% ============================================================
%% 辅助函数
%% ============================================================
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); if isempty(data), return; end
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w = prctile(data, 97.5);
    
    % 绘制所有数据点，而不仅仅是异常值，以增加散点数量
    if length(data) > 5000
        data_plot = data(randperm(length(data), 5000));
    else
        data_plot = data;
    end
    x_jit = center + (rand(size(data_plot)) - 0.5) * w_scatter;
    scatter(ax, x_jit, data_plot, sz_scatter, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
    
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    x_L = center - w_box/2; x_R = center + w_box/2;
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

function add_legend(ax, BoxBgColor, ScatterSize, OutlierText)
    h_box  = patch(ax, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
    h_med  = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2);
    h_mean = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    h_out  = scatter(ax, NaN, NaN, ScatterSize, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.6);
    legend(ax, [h_box, h_med, h_mean, h_out], {'IQR Box (25%-75%)', 'Median', 'Mean', OutlierText}, ...
           'Location', 'northeast', 'FontSize', 10, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');
end