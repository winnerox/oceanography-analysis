%% Analysis_MultiModel_Significance_Std.m
% =========================================================================
% 功能：三套数据 (EN4, IAP, Ishii) 标准态趋势与振幅的显著性检验
% 范围：1-8 阶
% 输出：MME Trend/Amp Maps (8 tiles), Boxplots
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
OutputDir = 'D:\work\Figures\StdRef';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Models = {'EN4', 'IAP', 'Ishii'};
Files_T = {'EN4_TSLA_Terms_1to8_StdRef.mat', ...
           'IAP_TSLA_Terms_1to8_StdRef.mat', ...
           'Ishii_TSLA_Terms_1to8_StdRef.mat'};
       
% Files_High removed (Consistency with Avg script: Analyze T-term only if HighOrder missing)

MaxOrder = 8; 

%% 2. 加载数据并计算
fprintf('>> [1/5] 加载三套数据并计算 Trend/Amp (Order 1-%d)...\n', MaxOrder);

Ref_Lon = []; Ref_Lat = [];
[Lon_Grid, Lat_Grid] = deal([]);

% 参考网格 (EN4)
d1 = load(fullfile(DataDir, Files_T{1}));
if isfield(d1, 'Lon'), Ref_Lon = d1.Lon; else, Ref_Lon = d1.lon; end
if isfield(d1, 'Lat'), Ref_Lat = d1.Lat; else, Ref_Lat = d1.lat; end
if max(Ref_Lon) > 180, Ref_Lon(Ref_Lon>180) = Ref_Lon(Ref_Lon>180)-360; end
[Ref_Lon, s_idx] = sort(Ref_Lon);
[Lon_Grid, Lat_Grid] = meshgrid(Ref_Lon, Ref_Lat);
[Ny, Nx] = size(Lon_Grid);

% 存储: [Ny, Nx, 3, 8]
All_Trend = zeros(Ny, Nx, 3, MaxOrder);
All_Amp   = zeros(Ny, Nx, 3, MaxOrder);

for m = 1:3
    fprintf('   Processing Model %d: %s ...\n', m, Models{m});
    
    dT = load(fullfile(DataDir, Files_T{m}));
    if isfield(dT, 'Lon'), cur_Lon = dT.Lon; else, cur_Lon = dT.lon; end
    if isfield(dT, 'Lat'), cur_Lat = dT.Lat; else, cur_Lat = dT.lat; end
    if max(cur_Lon) > 180, cur_Lon(cur_Lon>180) = cur_Lon(cur_Lon>180)-360; end
    [cur_Lon, sort_idx] = sort(cur_Lon);
    
    % HighOrder S/Cross term loading removed. Analysis is purely based on T terms.
    
    if isfield(dT, 'time_vec'), t_vec = dT.time_vec(:);
    elseif isfield(dT, 'time_axis'), t_vec = dT.time_axis(:);
    else, t_vec = dT.Time_Axis(:); 
    end
    
    NeedInterp = false;
    if length(cur_Lon) ~= length(Ref_Lon) || length(cur_Lat) ~= length(Ref_Lat) || ...
       any(abs(cur_Lon - Ref_Lon) > 1e-4) || any(abs(cur_Lat - Ref_Lat) > 1e-4)
        NeedInterp = true;
        [X_cur, Y_cur] = meshgrid(cur_Lon, cur_Lat);
    end
    
    for k = 1:MaxOrder
        % Term T
        Term_T = squeeze(dT.TSLA_AllOrders(sort_idx, :, :, k));
        
        % Term S/Cross ignored
        Total = Term_T;
        
        [tr_map, amp_map] = calc_trend_amp(Total, t_vec);
        
        if NeedInterp
            F_tr = griddedInterpolant(Y_cur, X_cur, tr_map', 'linear', 'none');
            F_am = griddedInterpolant(Y_cur, X_cur, amp_map', 'linear', 'none');
            tf = F_tr(Lat_Grid, Lon_Grid);
            af = F_am(Lat_Grid, Lon_Grid);
        else
            tf = tr_map'; af = amp_map';
        end
        All_Trend(:,:,m,k) = tf;
        All_Amp(:,:,m,k)   = af;
    end
end

%% 3. MME Stats
fprintf('>> [2/5] 计算 MME Stats ...\n');
Ens_Trend_Mean = squeeze(mean(All_Trend, 3)); 
Ens_Trend_Std  = squeeze(std(All_Trend, 0, 3));
Eq_SNR = abs(Ens_Trend_Mean) ./ Ens_Trend_Std; % SNR

Ens_Amp_Mean = squeeze(mean(All_Amp, 3));

%% 4. MME Trend Map (8 tiles)
fprintf('>> [3/5] 绘制 MME Trend Map ...\n');
figure('Position', [50, 50, 1200, 1400], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Diverging Cmap
MyBlue=[0.05 0.25 0.6]; MyLightBlue=[0.6 0.8 1]; MyWhite=[1 1 1];
MyLightRed=[1 0.7 0.6]; MyRed=[0.7 0.05 0.05];
ColorLocs=[0, 0.25, 0.5, 0.75, 1];
ColorPoints=[MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx=linspace(0,1,256);
MyCmap=[interp1(ColorLocs,ColorPoints(:,1),xx,'pchip')', interp1(ColorLocs,ColorPoints(:,2),xx,'pchip')', interp1(ColorLocs,ColorPoints(:,3),xx,'pchip')'];

for k = 1:MaxOrder
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Ens_Trend_Mean(:,:,k);
    
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % Dots (SNR > 1.69 for 90% CI with N=3)
    % t_crit(0.95, df=2) = 2.92. Threshold = 2.92 / sqrt(3) = 1.69
    mask_sig = Eq_SNR(:,:,k) > 1.69;
    skip = 3;
    [r, c] = find(mask_sig(1:skip:end, 1:skip:end));
    if ~isempty(r)
        ls = Lon_Grid(1:skip:end, 1:skip:end);
        las = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(ls(mask_sig(1:skip:end, 1:skip:end)), las(mask_sig(1:skip:end, 1:skip:end)), ...
               '.', 'color', [.2 .2 .2], 'markersize', 2);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
    colormap(gca, MyCmap);
    
    cl = prctile(abs(data(:)), 99.8); if cl<1e-5, cl=0.1; end
    caxis([-cl, cl]);
    
    title(sprintf('第 %d 阶 (Order %d)', k, k), 'FontSize', 12, 'FontWeight', 'bold');
    if mod(k,2)==0, cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8); end
end
sgtitle('\fontname{Microsoft YaHei}标准态 MME 趋势分布 (1-8阶, 点示 90% CI Significant)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Map_Trend_MME_Std.png'));

%% 5. MME Amp Map (8 tiles)
fprintf('>> [4/5] 绘制 MME Amp Map ...\n');
figure('Position', [100, 100, 1200, 1400], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Amp Cmap
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx_s=[0 0.1 0.4 0.8 1]; yy_s=linspace(0,1,256);
MySeqCmap=[interp1(xx_s,MyRedSeq(:,1),yy_s)', interp1(xx_s,MyRedSeq(:,2),yy_s)', interp1(xx_s,MyRedSeq(:,3),yy_s)'];

for k = 1:MaxOrder
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Ens_Amp_Mean(:,:,k);
    
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
    colormap(gca, MySeqCmap);
    
    cl = prctile(data(:), 99.8); if cl<1e-5, cl=0.1; end
    caxis([0, cl]);
    
    title(sprintf('第 %d 阶 (Order %d)', k, k), 'FontSize', 12, 'FontWeight', 'bold');
    if mod(k,2)==0, cb = colorbar; title(cb, 'mm', 'FontSize', 8); end
end
sgtitle('\fontname{Microsoft YaHei}标准态 MME 振幅分布 (1-8阶)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Map_Amp_MME_Std.png'));

%% 6. Boxplots
% Trend Boxplot (Figure 3)
% Flatten data for boxplot: [N_total, MaxOrder]
Trend_Flat = reshape(Ens_Trend_Mean, [], MaxOrder);
Amp_Flat   = reshape(Ens_Amp_Mean,   [], MaxOrder);

fprintf('>> [5.1/5] 绘制 MME Trend Boxplot (3-Axis Compact)...\n');
figure('Position', [300, 300, 1200, 550], 'Color', 'w', 'Name', 'MME Trend Compact');
MainPos = [0.07, 0.14, 0.80, 0.80];

% Parameters
C1 = [0.85 0.325 0.098]; % Red
C2 = [0 0.447 0.741];    % Blue
C3 = [0.466 0.674 0.188];% Green
BoxBgColor = [0.95 0.95 1];
BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.4; ScatterSize = 12;

% Data
Data_Trend = Trend_Flat;

% Layout
% AX1: 1-3
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:min(3, MaxOrder)
    draw_jitter_unit(ax1, k, Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = max(abs(prctile(Data_Trend(:,1:min(3,MaxOrder)), [0.2 99.8], 'all'))); 
if isnan(max_v1)||max_v1==0, max_v1=1; end
Limit1 = max_v1 * 1.1; ylim(ax1, [-Limit1, Limit1]);
set(ax1, 'XLim', [0.4, MaxOrder+1], 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% AX2: 4-5
ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 4:min(5, MaxOrder)
    draw_jitter_unit(ax2, k, Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
end
max_v2 = max(abs(prctile(Data_Trend(:,4:min(5,MaxOrder)), [0.2 99.8], 'all')));
if isnan(max_v2)||max_v2==0, max_v2=1; end
Limit2 = max_v2 * 1.2; ylim(ax2, [-Limit2, Limit2]);
set(ax2, 'Color', 'none', 'XLim', [0.4, MaxOrder+1], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% AX3: 6+
Offset_Ratio = 0.08; 
Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
ax3 = axes('Position', Pos3); hold(ax3, 'on');
Scale_Factor = Pos3(3) / MainPos(3);
XMin = 0.4;
XRange = (MaxOrder+1 - XMin) * Scale_Factor;
Real_X_Limit = [XMin, XMin + XRange];

for k = 6:MaxOrder
    draw_jitter_unit(ax3, k, Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
end
max_v3 = max(abs(prctile(Data_Trend(:,6:MaxOrder), [0.2 99.8], 'all')));
if isnan(max_v3)||max_v3==0, max_v3=1; end
Limit3 = max_v3 * 1.2; ylim(ax3, [-Limit3, Limit3]);
set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% Decoration
if ~isempty(ax1), yline(ax1, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.2); end
set(ax1, 'XTick', 1:MaxOrder, 'XTickLabel', arrayfun(@(x) sprintf('%d',x), 1:MaxOrder, 'UniformOutput', false));
xlabel(ax1, 'Order', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, 'Standard State MME Trend Distribution (Orders 1-8)', 'FontSize', 16, 'FontWeight', 'bold');
legend(ax1, {'Mean'}, 'Location', 'northeast', 'Box', 'off'); % simplified legend

saveas(gcf, fullfile(OutputDir, 'Stats_Trend_Boxplot_MME_Std.png'));


% Amp Boxplot (Figure 4)
fprintf('>> [5.2/5] 绘制 MME Amp Boxplot (3-Axis Compact)...\n');
figure('Position', [350, 350, 1200, 550], 'Color', 'w', 'Name', 'MME Amp Compact');
Data_Amp = Amp_Flat;

% AX1
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:min(3, MaxOrder)
    draw_jitter_unit(ax1, k, Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = prctile(Data_Amp(:,1:min(3,MaxOrder)), 99.8, 'all'); 
if isnan(max_v1)||max_v1==0, max_v1=1; end
ylim(ax1, [0, max_v1*1.1]);
set(ax1, 'XLim', [0.4, MaxOrder+1], 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% AX2
ax2 = axes('Position', MainPos); hold(ax2, 'on');
for k = 4:min(5, MaxOrder)
    draw_jitter_unit(ax2, k, Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
end
max_v2 = prctile(Data_Amp(:,4:min(5,MaxOrder)), 99.8, 'all'); 
if isnan(max_v2)||max_v2==0, max_v2=1; end
ylim(ax2, [0, max_v2*1.2]);
set(ax2, 'Color', 'none', 'XLim', [0.4, MaxOrder+1], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax2, 'Secondary Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% AX3
ax3 = axes('Position', Pos3); hold(ax3, 'on');
Scale_Factor = Pos3(3) / MainPos(3);
XMin = 0.4;
XRange = (MaxOrder+1 - XMin) * Scale_Factor;
Real_X_Limit = [XMin, XMin + XRange];

for k = 6:MaxOrder
    draw_jitter_unit(ax3, k, Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
end
max_v3 = prctile(Data_Amp(:,6:MaxOrder), 99.8, 'all'); 
if isnan(max_v3)||max_v3==0, max_v3=1; end
ylim(ax3, [0, max_v3*1.2]);
set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax3, 'High-Order Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

set(ax1, 'XTick', 1:MaxOrder, 'XTickLabel', arrayfun(@(x) sprintf('%d',x), 1:MaxOrder, 'UniformOutput', false));
xlabel(ax1, 'Order', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, 'Standard State MME Amplitude Distribution (Orders 1-8)', 'FontSize', 16, 'FontWeight', 'bold');
legend(ax1, {'Mean'}, 'Location', 'northeast', 'Box', 'off');

saveas(gcf, fullfile(OutputDir, 'Stats_Amp_Boxplot_MME_Std.png'));

fprintf('>> Done!\n');
%% Internal
function [t_map, a_map] = calc_trend_amp(data_3d, t_vec)
    [nx, ny, nt] = size(data_3d);
    t_map = nan(nx, ny); a_map = nan(nx, ny);
    for i=1:nx
        for j=1:ny
            y = squeeze(data_3d(i,j,:));
            idx = ~isnan(y);
            if sum(idx) < nt*0.6, continue; end
            p = polyfit(t_vec(idx), y(idx), 1);
            t_map(i,j) = p(1);
            resid = y(idx) - polyval(p, t_vec(idx));
            a_map(i,j) = std(resid);
        end
    end
end

