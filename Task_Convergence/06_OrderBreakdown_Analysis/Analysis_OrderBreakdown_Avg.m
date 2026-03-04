%% Analysis_OrderBreakdown_Avg.m
% =========================================================================
% 功能：平均态各阶贡献分解 (1阶, 2阶, 3阶)
% 展示：各阶的趋势和振幅 空间分布 + 箱线图
% 风格：与 SingleModel_Plots 保持一致 (Robinson投影 + Jitter箱线图)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
OutputDir = 'D:\work\Figures\OrderBreakdown';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> [1/4] 加载数据...\n');

% 纯T项
f_T = fullfile(DataDir, 'EN4_TSLA_Terms_1to8_Average.mat');
if ~exist(f_T, 'file'), error('找不到纯T项数据'); end
D_T = load(f_T);

% 经纬度/时间轴修正
if isfield(D_T, 'Lon'), Lon = D_T.Lon; elseif isfield(D_T, 'longitude'), Lon = D_T.longitude; end
if isfield(D_T, 'Lat'), Lat = D_T.Lat; elseif isfield(D_T, 'latitude'), Lat = D_T.latitude; end
if max(Lon) > 180, Lon(Lon > 180) = Lon(Lon > 180) - 360; end
[Lon, sort_idx] = sort(Lon);

% 稳健的时间变量读取
if isfield(D_T, 'time_vec'), time_vec = D_T.time_vec;
elseif isfield(D_T, 'time_axis'), time_vec = D_T.time_axis;
elseif isfield(D_T, 'Time_Axis'), time_vec = D_T.Time_Axis;
else, error('无法找到时间变量'); end
t_vec = time_vec(:);

TSLA_AllOrders = D_T.TSLA_AllOrders(sort_idx, :, :, :);
T1 = squeeze(TSLA_AllOrders(:,:,:,1));
T2 = squeeze(TSLA_AllOrders(:,:,:,2));
T3 = squeeze(TSLA_AllOrders(:,:,:,3));

%% 3. 计算统计量 (1-3阶)
fprintf('>> [2/4] 计算统计量 (1-3阶)...\n');
[Nx, Ny, Nt] = size(T1);
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

% 计算函数 (Robust Loop)
calc_stats = @(data) compute_stats_robust(data, t_vec);

fprintf('   计算 Order 1...\n'); [Trend1, Amp1, ~] = calc_stats(T1);
fprintf('   计算 Order 2...\n'); [Trend2, Amp2, ~] = calc_stats(T2);
fprintf('   计算 Order 3...\n'); [Trend3, Amp3, ~] = calc_stats(T3);

% 智能单位修正 (对齐 SingleModel_Plots 逻辑)
% 如果 Trend > 500 (mm/yr)，说明可能是被错误乘了 1000 (即已经是 mm 但又乘了)
% 或者如果数据源是 m，那 calc_stats 里没乘 1000，这里需要乘。
% 现在的 calc_stats 不乘 1000。
% 所以如果 max(Trend) < 0.1 (mm/yr 很小)，可能是 m，自动转换为 mm。
% Unit Heuristic removed (Risky). Assuming input data is already correctly scaled.

%% 4. 绘图 - 空间分布 (Trends)
fprintf('>> [3/4] 绘制 Trend 空间分布...\n');

Data_Map = {Trend1, Trend2, Trend3};
Titles = {'1阶趋势', '2阶趋势', '3阶趋势'};

% Colormap (Diverging for Trend)
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];

figure('Position', [100, 100, 1200, 400], 'Color', 'w');
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:3
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid, Lat_Grid, Data_Map{k}); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8);
    colormap(gca, MyDivergentCmap);
    
    % Auto-scale 99%
    clim = prctile(abs(Data_Map{k}(:)), 99);
    if isnan(clim) || clim == 0, clim = 0.5; end
    caxis([-clim, clim]);
    
    title(Titles{k}, 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8);
end
sgtitle('\fontname{Microsoft YaHei}各阶趋势分布 (1-3阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Trend_Order123.png'));

%% 5. 绘图 - 空间分布 (Amplitudes)
fprintf('>> [4/4] 绘制 Amplitude 空间分布...\n');

Data_Amp = {Amp1, Amp2, Amp3};
Titles_Amp = {'1阶振幅', '2阶振幅', '3阶振幅'};

% Colormap (Sequential for Amp)
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx_s = [0 0.1 0.4 0.8 1]; yy_s = linspace(0,1,256);
MySeqCmap = [interp1(xx_s,MyRedSeq(:,1),yy_s)', interp1(xx_s,MyRedSeq(:,2),yy_s)', interp1(xx_s,MyRedSeq(:,3),yy_s)'];

figure('Position', [100, 600, 1200, 400], 'Color', 'w');
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:3
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid, Lat_Grid, Data_Amp{k}); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8);
    colormap(gca, MySeqCmap);
    
    % Auto-scale 99%
    clim = prctile(Data_Amp{k}(:), 99);
    if isnan(clim) || clim == 0, clim = 1; end
    caxis([0, clim]);
    
    title(Titles_Amp{k}, 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; title(cb, 'mm', 'FontSize', 8);
end
sgtitle('\fontname{Microsoft YaHei}各阶振幅分布 (1-3阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Amp_Order123.png'));

%% 6. 箱线图 - Trend (1-3阶对比)
fprintf('>> [5/6] 绘制 Trend 箱线图...\n');

Data_Box_T = {Trend1, Trend2, Trend3};
Labels_T = {'1阶', '2阶', '3阶'};
Colors_T = {[0.85 0.325 0.098], [0 0.447 0.741], [0.466 0.674 0.188]};

figure('Position', [100, 100, 500, 400], 'Color', 'w');
hold on;
BoxWidth = 0.5; ScatterWidth = 0.3;

    data = Data_Box_T{k}(:);
    draw_jitter_unit(gca, k, data, BoxWidth, ScatterWidth, 10, Colors_T{k}, [0.95 0.95 1], 0.4);
end
% Manual Legend (Keep existing)
h_box = patch(NaN, NaN, [0.95 0.95 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
h_med = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
legend([h_box, h_med, h_mean], {'IQR (25%-75%)', 'Median', 'Mean'}, ...
       'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
set(gca, 'XTick', 1:3, 'XTickLabel', Labels_T);
ylabel('Trend (mm/yr)', 'FontSize', 12, 'FontWeight', 'bold');
title('\fontname{Microsoft YaHei}各阶趋势对比 (1-3阶)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
saveas(gcf, fullfile(OutputDir, 'Boxplot_Trend_Order123.png'));

%% 7. 箱线图 - Amplitude (1-3阶对比)
fprintf('>> [6/6] 绘制 Amplitude 箱线图...\n');

Data_Box_A = {Amp1, Amp2, Amp3};

figure('Position', [650, 100, 500, 400], 'Color', 'w');
hold on;

for k = 1:3
    data = Data_Box_A{k}(:); 
    draw_jitter_unit(gca, k, data, BoxWidth, ScatterWidth, 10, Colors_T{k}, [0.95 0.95 1], 0.4);
end
h_box = patch(NaN, NaN, [0.95 0.95 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
h_med = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
legend([h_box, h_med, h_mean], {'IQR (25%-75%)', 'Median', 'Mean'}, ...
       'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
set(gca, 'XTick', 1:3, 'XTickLabel', Labels_T);
ylabel('Amplitude (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title('\fontname{Microsoft YaHei}各阶振幅对比 (1-3阶)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
saveas(gcf, fullfile(OutputDir, 'Boxplot_Amp_Order123.png'));

fprintf('\n>> 分析完成！图像已保存至 %s\n', OutputDir);

%% Helper Function
function [trend, amp, sig] = compute_stats_robust(data_3d, t_vec)
    % Robust loop-based calculation
    [Nx, Ny, Nt] = size(data_3d);
    
    % Ensure column vector
    t_vec = t_vec(:); 
    
    trend = nan(Nx, Ny); 
    amp = nan(Nx, Ny); 
    sig = nan(Nx, Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            y = squeeze(data_3d(i, j, :));
            mask = ~isnan(y);
            if sum(mask) < Nt * 0.7, continue; end
            
            x_sub = t_vec(mask);
            y_sub = y(mask);
            
            % Trend
            p = polyfit(x_sub, y_sub, 1);
            trend(i, j) = p(1);
            
            % Amp
            y_fit = polyval(p, x_sub);
            resid = y_sub - y_fit;
            amp(i, j) = std(resid);
            
            % Sig (90% Confidence Interval -> alpha = 0.10)
            [~, ~, ~, ~, stats] = regress(y_sub, [ones(length(x_sub),1), x_sub]);
            if stats(3) < 0.10, sig(i, j) = 1; else, sig(i, j) = 0; end
        end
    end
    
    % Transpose to (Lat, Lon) for m_pcolor
    trend = trend';
    amp = amp';
    sig = sig';
end
