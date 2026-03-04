%% Analysis_Convergence.m
% =========================================================================
% 功能：平均态收敛性分析 (Average State Convergence Analysis)
% 比较：
%   1. 精确解 (Exact): 全非线性计算，无近似
%   2. 3阶近似 (Approximation): T(1-3) + S(1-3) + Cross(2-3)
%   3. 残差 (Residual): Exact - Approx
% 风格：与 SingleModel_Plots 保持一致 (Robinson投影 + Jitter箱线图)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
TermsDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\Convergence';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> [1/5] 加载数据...\n');

% 2.1 精确解
f_exact = fullfile(TermsDir, 'EN4_Formula11_Exact_Avg.mat');
if ~exist(f_exact, 'file'), error('找不到精确解数据: %s', f_exact); end
D1 = load(f_exact);
TSLA_Exact = D1.TSLA_Exact_Avg; 

% 2.1 精确解
TSLA_Exact_Avg = []; load(fullfile(TermsDir, 'EN4_Formula11_Exact_Avg.mat')); 
% 假设变量名为 TSLA_Exact_Avg, 若不同需调整
if exist('TSLA_Exact_Avg','var'), TSLA_Exact = TSLA_Exact_Avg; end

% 稳健的时间变量读取
if isfield(D1, 'time_vec')
    time_vec = D1.time_vec;
elseif isfield(D1, 'time_axis')
    time_vec = D1.time_axis;
elseif isfield(D1, 'Time_Axis')
    time_vec = D1.Time_Axis;
else
    error('无法找到时间变量 (time_vec / time_axis)');
end

Lon = D1.lon; Lat = D1.lat;
if max(Lon) > 180, Lon(Lon > 180) = Lon(Lon > 180) - 360; end
[Lon, sort_idx] = sort(Lon);
TSLA_Exact = TSLA_Exact(sort_idx, :, :);

% 2.2 纯T项 (1-3阶)
f_T = fullfile(DataDir, 'EN4_TSLA_Terms_1to8_Average.mat');
if ~exist(f_T, 'file'), error('找不到纯T项数据'); end
D2 = load(f_T);
TSLA_AllOrders = D2.TSLA_AllOrders(sort_idx, :, :, :);
T1 = squeeze(TSLA_AllOrders(:,:,:,1));
T2 = squeeze(TSLA_AllOrders(:,:,:,2));
T3 = squeeze(TSLA_AllOrders(:,:,:,3));

% 2.3 盐度项与混合项
f_S = fullfile(TermsDir, 'EN4_Terms_S_Cross_Avg_1to3.mat');
if ~exist(f_S, 'file'), error('找不到S项及混合项数据'); end
D3 = load(f_S);
S1 = D3.SSLA_1(sort_idx, :, :); 
S2 = D3.SSLA_2(sort_idx, :, :); 
S3 = D3.SSLA_3(sort_idx, :, :);
C_ST = D3.Cross_ST(sort_idx, :, :); 
C_SST = D3.Cross_SST(sort_idx, :, :); 
C_STT = D3.Cross_STT(sort_idx, :, :);

%% 3. 计算近似解与残差
fprintf('>> [2/5] 计算近似解与残差...\n');

% 3阶近似 = T(1-3) + S(1-3) + Cross(ST, SST, STT)
TSLA_Approx = (T1 + S1) + (T2 + S2 + C_ST) + (T3 + S3 + C_SST + C_STT);

Residual = TSLA_Exact - TSLA_Approx;

%% 4. 统计趋势与振幅
fprintf('>> [3/5] 计算趋势与振幅...\n');
[Nx, Ny, Nt] = size(TSLA_Exact);
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

% 数据已在 Calc 脚本中转换为 mm，无需再做单位转换

% 计算 Function Handle
calc_stats = @(data) compute_stats_robust(data, time_vec);

fprintf('   计算 Residual 统计量...\n');
[Trend_Res, Amp_Res, ~] = calc_stats(Residual);

fprintf('   计算 Exact 统计量...\n');
[Trend_Exact, ~, ~] = calc_stats(TSLA_Exact);

fprintf('   计算 Approx 统计量...\n');
[Trend_Approx, ~, ~] = calc_stats(TSLA_Approx);


%% 5. 绘图 - 残差分布图 (Robinson + Diverging)
fprintf('>> [4/5] 绘制残差空间分布...\n');

figure('Position', [100, 100, 1000, 450], 'Color', 'w');
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Colormaps
% Diverging for Trend
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];
% Sequential for Amp
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx_s = [0 0.1 0.4 0.8 1]; yy_s = linspace(0,1,256);
MySeqCmap = [interp1(xx_s,MyRedSeq(:,1),yy_s)', interp1(xx_s,MyRedSeq(:,2),yy_s)', interp1(xx_s,MyRedSeq(:,3),yy_s)'];

% --- Plot Trend Residual ---
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, Trend_Res); shading interp;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
       'linewidth', 0.5, 'fontsize', 8);
colormap(gca, MyDivergentCmap);
% Smart auto-scale: use 95th percentile, round to clean number
clim_val = prctile(abs(Trend_Res(:)), 95);
if isnan(clim_val) || clim_val == 0, clim_val = 0.5; end
clim_val = ceil(clim_val * 1000) / 1000; % round up to 3 decimal places
caxis([-clim_val, clim_val]);

title('平均态残差趋势 (Residual Trend)', 'FontSize', 12, 'FontWeight', 'bold');
cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8);

% --- Plot Amp Residual ---
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, Amp_Res); shading interp;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
       'linewidth', 0.5, 'fontsize', 8);
colormap(gca, MySeqCmap);
% Smart auto-scale for Amp
amp_95 = prctile(Amp_Res(:), 95);
if isnan(amp_95) || amp_95 == 0, amp_95 = 1; end
amp_95 = ceil(amp_95 * 100) / 100; % round up to 2 decimal places
caxis([0, amp_95]);

title('平均态残差振幅 (Residual Amplitude)', 'FontSize', 12, 'FontWeight', 'bold');
cb = colorbar; title(cb, 'mm', 'FontSize', 8);

sgtitle('\fontname{Microsoft YaHei}平均态残差分布 (Exact - 3rd Order Approx)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Residual_Map.png'));

%% 6. 绘图 - 收敛性验证 (Jitter Boxplot)
fprintf('>> [5/5] 绘制收敛性检验箱线图...\n');

Data_Box = {Trend_Exact, Trend_Approx, Trend_Res};
Labels = {'Exact', 'Approx (3rd)', 'Residual'};
Colors = {[0.2 0.2 0.2], [0 0.447 0.741], [0.85 0.325 0.098]};

figure('Position', [100, 100, 600, 450], 'Color', 'w');
hold on;

BoxWidth = 0.5;
ScatterWidth = 0.3;

for k = 1:3
    data = Data_Box{k}(:);
    data = data(~isnan(data));
    if isempty(data), continue; end
    
    center = k;
    col = Colors{k};
    
    % Stats
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
    
    % Whiskers
    plot([center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot([center, center], [q3, up_w], 'k-', 'LineWidth', 1.2);
    plot([center-BoxWidth*0.2, center+BoxWidth*0.2], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot([center-BoxWidth*0.2, center+BoxWidth*0.2], [up_w, up_w], 'k-', 'LineWidth', 1.2);
    
    % Mean/Median
    plot([center-BoxWidth/2, center+BoxWidth/2], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

% Manual Legend
h_box = patch(NaN, NaN, [0.95 0.95 1], 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
h_med = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_sca = scatter(NaN, NaN, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.4);

legend([h_box, h_med, h_mean, h_sca], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 10, 'Box', 'on');

set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
ylabel('趋势 (Trend) [mm/yr]', 'FontSize', 12, 'FontWeight', 'bold');
title('平均态收敛性检验: 精确解 vs 近似解', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([-5 5]); 

saveas(gcf, fullfile(OutputDir, 'Convergence_Boxplot.png'));

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
            % At least 70% valid data required
            if sum(mask) < Nt * 0.7, continue; end
            
            x_sub = t_vec(mask);
            y_sub = y(mask);
            
            % Trend (Polyfit)
            p = polyfit(x_sub, y_sub, 1);
            trend(i, j) = p(1);
            
            % Amplitude (Std of Detrended)
            y_fit = polyval(p, x_sub);
            resid = y_sub - y_fit;
            amp(i, j) = std(resid);
            
            % Significance (Regress p-value)
            % [b, bint, r, rint, stats] = regress(y, X)
            [~, ~, ~, ~, stats] = regress(y_sub, [ones(length(x_sub),1), x_sub]);
            % stats(3) is p-value
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
