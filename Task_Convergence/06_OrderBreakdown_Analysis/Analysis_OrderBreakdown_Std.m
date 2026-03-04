%% Analysis_OrderBreakdown_Std.m
% =========================================================================
% 功能：标准态 (Std Ref) 1-8 阶收敛性分解分析
% 内容：
%   1. 计算每一阶的总贡献 (Term_k = T_k + S_k + Cross_k)
%   2. 绘制 1-8 阶的 Trend 空间分布 (4x2 Tiled Layout)
%   3. 绘制 1-8 阶的 Total Amplitude 空间分布 (4x2 Tiled Layout)
%   4. 绘制各阶 Trend 和 Amplitude 的箱线图
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
TermsDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\StdRef';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> [1/4] 加载数据 (1-8阶)... \n');

% 2.1 Pure T (1-8阶)
f_T = fullfile(DataDir, 'EN4_TSLA_Terms_1to8_StdRef.mat');
if ~exist(f_T, 'file'), error('找不到纯T项数据: %s', f_T); end
D_T = load(f_T);
% 健壮读取 Lon/Lat (兼容不同命名)
if isfield(D_T, 'Lon'), Lon = D_T.Lon;
elseif isfield(D_T, 'lon'), Lon = D_T.lon;
elseif isfield(D_T, 'longitude'), Lon = D_T.longitude;
else, error('T项数据中找不到经度变量'); end

if isfield(D_T, 'Lat'), Lat = D_T.Lat;
elseif isfield(D_T, 'lat'), Lat = D_T.lat;
elseif isfield(D_T, 'latitude'), Lat = D_T.latitude;
else, error('T项数据中找不到纬度变量'); end

if max(Lon) > 180, Lon(Lon > 180) = Lon(Lon > 180) - 360; end
[Lon, sort_idx] = sort(Lon);

% Load T Terms [Nx, Ny, Time, 8]
TSLA_T = D_T.TSLA_AllOrders(sort_idx, :, :, 1:8);

% 2.2 Pure S + Cross (1-8阶)
f_S = fullfile(TermsDir, 'EN4_Std_Terms_HighOrder_1to8.mat');
if ~exist(f_S, 'file'), error('找不到S项及混合项数据: %s', f_S); end
D_S = load(f_S);

% 时间变量
if isfield(D_T, 'time_vec'), t_vec = D_T.time_vec;
elseif isfield(D_T, 'time_axis'), t_vec = D_T.time_axis;
elseif isfield(D_T, 'Time_Axis'), t_vec = D_T.Time_Axis;
elseif isfield(D_S, 'time_vec'), t_vec = D_S.time_vec;
else, error('无法找到时间变量'); end
t_vec = t_vec(:);

%% 3. 逐阶计算统计量
fprintf('>> [2/4] 计算各阶 Trend 和 Amplitude...\n');

[Nx, Ny, ~] = size(TSLA_T, 1:3);
Lat_Grid = []; Lon_Grid = [];
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

Trend_Maps = cell(1, 8);
Amp_Maps = cell(1, 8);

calc_stats = @(data) compute_stats_vectorized(data, t_vec);

for k = 1:8
    fprintf('   Processing Order %d ... ', k);
    
    % T Term
    Term_k = squeeze(TSLA_T(:,:,:,k));
    
    % S Term
    fname_s = sprintf('SSLA_%d', k);
    if isfield(D_S, fname_s)
        tmp_s = D_S.(fname_s)(sort_idx, :, :);
        Term_k = Term_k + tmp_s;
    end
    
    % Cross Term (k>=2)
    if k >= 2
        fname_c = sprintf('Cross_%d', k);
        if isfield(D_S, fname_c)
            tmp_c = D_S.(fname_c)(sort_idx, :, :);
            Term_k = Term_k + tmp_c;
        end
    end
    
    % Compute Stats
    [tr, am, ~] = calc_stats(Term_k);
    Trend_Maps{k} = tr;
    Amp_Maps{k} = am;
    
    fprintf('Done.\n');
end

%% 4. 绘图 - Trend Maps (4x2)
fprintf('>> [3/4] 绘制 Trend 分布 (Order 1-8)...\n');

% Colormap
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];

figure('Position', [50, 50, 1600, 800], 'Color', 'w');
t = tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:8
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid, Lat_Grid, Trend_Maps{k}); shading interp;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 7, 'xticklabels', [], 'yticklabels', []);
    colormap(gca, MyDivergentCmap);
    
    % 智能设置 Colorbar 范围 (98% 分位数，避免极端值)
    vals = Trend_Maps{k}(:);
    clim_val = prctile(abs(vals), 98);
    if isnan(clim_val) || clim_val < 0.01, clim_val = 0.1; end
    clim_val = ceil(clim_val * 10) / 10;
    caxis([-clim_val, clim_val]);
    
    title(sprintf('第 %d 阶趋势', k), 'FontSize', 11, 'FontWeight', 'bold');
    cb = colorbar; 
    cb.Label.String = 'mm/yr';
    cb.FontSize = 7;
end
sgtitle('\fontname{Microsoft YaHei}标准态各阶趋势分布 (1-8阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Trend_Maps_Std_Order1to8.png'));

%% 5. 绘图 - Amp Maps (4x2)
fprintf('>> [4/4] 绘制 Amplitude 分布 (Order 1-8)...\n');

MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx_s = [0 0.1 0.4 0.8 1]; yy_s = linspace(0,1,256);
MySeqCmap = [interp1(xx_s,MyRedSeq(:,1),yy_s)', interp1(xx_s,MyRedSeq(:,2),yy_s)', interp1(xx_s,MyRedSeq(:,3),yy_s)'];

figure('Position', [100, 100, 1600, 800], 'Color', 'w');
t = tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:8
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid, Lat_Grid, Amp_Maps{k}); shading interp;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 7, 'xticklabels', [], 'yticklabels', []);
    colormap(gca, MySeqCmap);
    
    % Smart Caxis
    vals = Amp_Maps{k}(:);
    clim_val = prctile(vals, 98);
    if isnan(clim_val) || clim_val < 0.01, clim_val = 0.5; end
    clim_val = ceil(clim_val * 10) / 10;
    caxis([0, clim_val]);
    
    title(sprintf('第 %d 阶振幅', k), 'FontSize', 11, 'FontWeight', 'bold');
    cb = colorbar; title(cb, 'mm', 'FontSize', 7);
end
sgtitle('\fontname{Microsoft YaHei}标准态各阶振幅分布 (1-8阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Amp_Maps_Std_Order1to8.png'));

%% 6. 箱线图 (Trend & Amp Summary)
figure('Position', [100, 100, 1000, 450], 'Color', 'w');
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Trend Boxplot
nexttile; hold on;
BoxWidth = 0.6; 
Colors = jet(8);
for k = 1:8
    data = Trend_Maps{k}(:); 
    draw_jitter_unit(gca, k, data, BoxWidth, 0.3, 10, Colors(k,:), [0.95 0.95 1], 0.6);
end
xlim([0.5, 8.5]); xticks(1:8); xlabel('Order'); ylabel('Trend (mm/yr)');
title('标准态各阶趋势对比', 'FontSize', 12); grid on;

% Amp Boxplot
nexttile; hold on;
for k = 1:8
    data = Amp_Maps{k}(:); 
    draw_jitter_unit(gca, k, data, BoxWidth, 0.3, 10, Colors(k,:), [0.95 0.95 1], 0.6);
end
xlim([0.5, 8.5]); xticks(1:8); xlabel('Order'); ylabel('Amplitude (mm)');
title('标准态各阶振幅对比', 'FontSize', 12); grid on;

sgtitle('\fontname{Microsoft YaHei}标准态 1-8 阶统计分布概览', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Boxplot_Summary_Std_Order1to8.png'));

fprintf('\n>> 分析完成！图像已保存至 %s\n', OutputDir);

%% Helper Function
function [trend, amp, sig] = compute_stats_vectorized(data_3d, t_vec)
    [Nx, Ny, Nt] = size(data_3d);
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
            
            p = polyfit(x_sub, y_sub, 1);
            trend(i, j) = p(1);
            
            y_fit = polyval(p, x_sub);
            resid = y_sub - y_fit;
            amp(i, j) = std(resid);
            
            [~, ~, ~, ~, stats] = regress(y_sub, [ones(length(x_sub),1), x_sub]);
            if stats(3) < 0.10 % 90% CI
                sig(i, j) = 1; else, sig(i, j) = 0; end
        end
    end
    trend = trend'; amp = amp'; sig = sig';
end
