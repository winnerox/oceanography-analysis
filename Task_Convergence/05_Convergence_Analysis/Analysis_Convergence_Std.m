%% Analysis_Convergence_Std.m
% =========================================================================
% 功能：标准态收敛性分析 (Std Ref State Convergence Analysis)
% 比较：
%   1. 精确解 (Exact): 全非线性计算，无近似
%   2. 8阶近似 (Approximation): T(1-8) + S(1-8) + Cross(2-8)
%   3. 残差 (Residual): Exact - Approx
% 风格：与 SingleModel_Plots 保持一致 (Robinson投影 + Jitter箱线图)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
TermsDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\StdRef';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> [1/5] 加载数据 (标准态 8阶)...\n');

% 2.1 精确解
f_exact = fullfile(TermsDir, 'EN4_Formula11_Exact_Std.mat');
if ~exist(f_exact, 'file'), error('找不到精确解数据: %s', f_exact); end
D1 = load(f_exact);
TSLA_Exact = D1.TSLA_Exact_Std; 
time_vec = D1.time_vec;
Lon = D1.lon; Lat = D1.lat;
if max(Lon) > 180, Lon(Lon > 180) = Lon(Lon > 180) - 360; end
[Lon, sort_idx] = sort(Lon);
TSLA_Exact = TSLA_Exact(sort_idx, :, :);

% 2.2 Pure T (1-8阶)
f_T = fullfile(DataDir, 'EN4_TSLA_Terms_1to8_StdRef.mat');
if ~exist(f_T, 'file'), error('找不到纯T项数据'); end
D2 = load(f_T);
TSLA_AllOrders = D2.TSLA_AllOrders(sort_idx, :, :, :);
% Sum T 1-8
Sum_T = sum(TSLA_AllOrders(:,:,:,1:8), 4);

% 2.3 Pure S + Cross (1-8阶)
f_S = fullfile(TermsDir, 'EN4_Std_Terms_HighOrder_1to8.mat');
if ~exist(f_S, 'file'), error('找不到S项及混合项数据'); end
D3 = load(f_S);

% 提取并对应排序
% Pure S
Sum_S = zeros(size(TSLA_Exact));
for k=1:8
    fname = sprintf('SSLA_%d', k);
    if isfield(D3, fname)
        tmp = D3.(fname)(sort_idx, :, :);
        Sum_S = Sum_S + tmp;
    end
end

% Cross Terms (只有 2-8 阶)
Sum_Cross = zeros(size(TSLA_Exact));
for k=2:8
    fname = sprintf('Cross_%d', k);
    if isfield(D3, fname)
        tmp = D3.(fname)(sort_idx, :, :); 
        Sum_Cross = Sum_Cross + tmp;
    end
end

%% 3. 计算近似解与残差
fprintf('>> [2/5] 计算近似解与残差...\n');

% 8阶近似
TSLA_Approx = Sum_T + Sum_S + Sum_Cross;
Residual = TSLA_Exact - TSLA_Approx;

%% 4. 统计趋势与振幅
fprintf('>> [3/5] 计算趋势与振幅...\n');
[Nx, Ny, Nt] = size(TSLA_Exact);
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

t_vec = time_vec(:);

% 向量化计算统计量
calc_stats = @(data) compute_stats_vectorized(data, t_vec);

fprintf('   计算 Residual 统计量...\n');
[Trend_Res, Amp_Res, ~] = calc_stats(Residual);

fprintf('   计算 Exact 统计量...\n');
[Trend_Exact, ~, ~] = calc_stats(TSLA_Exact);

fprintf('   计算 Approx 统计量...\n');
[Trend_Approx, ~, ~] = calc_stats(TSLA_Approx);

%% 5. 绘图 - 残差分布图 (Robinson + Diverging + Interp)
fprintf('>> [4/5] 绘制残差空间分布...\n');

% Colormaps
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx_s = [0 0.1 0.4 0.8 1]; yy_s = linspace(0,1,256);
MySeqCmap = [interp1(xx_s,MyRedSeq(:,1),yy_s)', interp1(xx_s,MyRedSeq(:,2),yy_s)', interp1(xx_s,MyRedSeq(:,3),yy_s)'];

figure('Position', [100, 100, 1000, 450], 'Color', 'w');
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot Trend Residual
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, Trend_Res); shading interp;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
colormap(gca, MyDivergentCmap);

% Smart Auto-scale for Trend Residual (95% round up)
% Smart Auto-scale for Trend Residual (95% round up)
clim_val = prctile(abs(Trend_Res(:)), 98);
fprintf('   [Debug] Trend Residual Range: Min=%.6f, Max=%.6f, Mean=%.6f\n', min(Trend_Res(:)), max(Trend_Res(:)), nanmean(Trend_Res(:)));
if isnan(clim_val) || clim_val == 0, clim_val = 1e-5; end % Valid small number
% clim_val = ceil(clim_val * 100) / 100;  % Remove rounding to see small values
caxis([-clim_val, clim_val]);

title('标准态残差趋势 (8阶)', 'FontSize', 12, 'FontWeight', 'bold');
cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8);

% Plot Amp Residual
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, Amp_Res); shading interp;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
colormap(gca, MySeqCmap);

% Smart Auto-scale for Amp Residual (95% round up)
% Smart Auto-scale for Amp Residual (95% round up)
amp_95 = prctile(Amp_Res(:), 98);
fprintf('   [Debug] Amp Residual Range: Min=%.6f, Max=%.6f, Mean=%.6f\n', min(Amp_Res(:)), max(Amp_Res(:)), nanmean(Amp_Res(:)));
if isnan(amp_95) || amp_95 == 0, amp_95 = 1e-4; end
% amp_95 = ceil(amp_95 * 10) / 10;
caxis([0, amp_95]);

title('标准态残差振幅', 'FontSize', 12, 'FontWeight', 'bold');
cb = colorbar; title(cb, 'mm', 'FontSize', 8);

sgtitle('\fontname{Microsoft YaHei}标准态 8阶收敛性残差 (Exact - 8th Approx)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Residual_Map_Std_Order8.png'));

%% 6. 绘图 - 收敛性验证 (Jitter Boxplot)
fprintf('>> [5/5] 绘制收敛性检验箱线图...\n');

Data_Box = {Trend_Exact, Trend_Approx, Trend_Res};
Labels = {'Exact (Std)', '8th Order Approx', 'Residual'};
Colors = {[0.2 0.2 0.2], [0 0.447 0.741], [0.85 0.325 0.098]};

figure('Position', [100, 100, 600, 450], 'Color', 'w');
hold on;
BoxWidth = 0.5; ScatterWidth = 0.3;

for k = 1:3
    data = Data_Box{k}(:); data = data(~isnan(data));
    if isempty(data), continue; end
    
    center = k; col = Colors{k};
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w = prctile(data, 97.5);
    
    % Jitter Outliers
    idx_draw = data < low_w | data > up_w;
    vals_draw = data(idx_draw);
    if length(vals_draw) > 1000, idx_sub=randperm(length(vals_draw),1000); vals_draw=vals_draw(idx_sub); end
    x_jit = center + (rand(size(vals_draw))-0.5)*ScatterWidth;
    scatter(x_jit, vals_draw, 10, col, 'filled', 'MarkerFaceAlpha', 0.4);
    
    % Box
    patch([center-BoxWidth/2, center+BoxWidth/2, center+BoxWidth/2, center-BoxWidth/2], ...
          [q1, q1, q3, q3], [0.95 0.95 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
    plot([center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot([center, center], [q3, up_w], 'k-', 'LineWidth', 1.2);
    plot([center-BoxWidth*0.2, center+BoxWidth*0.2], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot([center-BoxWidth*0.2, center+BoxWidth*0.2], [up_w, up_w], 'k-', 'LineWidth', 1.2);
    plot([center-BoxWidth/2, center+BoxWidth/2], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

h_box = patch(NaN, NaN, [0.95 0.95 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
h_med = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
legend([h_box, h_med, h_mean], {'IQR (25%-75%)', 'Median', 'Mean'}, ...
       'Location', 'northwest', 'FontSize', 10, 'Box', 'on');

set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
ylabel('趋势 (Trend) [mm/yr]', 'FontSize', 12, 'FontWeight', 'bold');
title('标准态收敛性检验', 'FontSize', 14, 'FontWeight', 'bold');
grid on; ylim([-5 5]); % Adjust y-limit as needed
saveas(gcf, fullfile(OutputDir, 'Convergence_Boxplot_Std_Order8.png'));

fprintf('\n>> 分析完成！图像已保存至 %s\n', OutputDir);

%% Helper Function
function [trend, amp, sig] = compute_stats_vectorized(data_3d, t_vec)
    % Robust loop-based calculation
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
                sig(i, j) = 1; 
            else
                sig(i, j) = 0;
            end
        end
    end
    
    % Transpose to (Lat, Lon) for m_pcolor
    trend = trend';
    amp = amp';
    sig = sig';
end
