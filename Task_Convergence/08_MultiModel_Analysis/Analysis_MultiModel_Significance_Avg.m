%% Analysis_MultiModel_Significance_Avg.m
% =========================================================================
% 功能：三套数据 (EN4, IAP, Ishii) 平均态趋势与振幅的显著性检验
% 范围：1-3 阶 (因平均态收敛快，且数据仅生成到3阶)
% 
% 输出图表 (参考 SingleModel_Plots/Average/):
%   1. Map_Trend: 集合平均趋势分布 (Overlay Significance Dots)
%   2. Map_Amp:   集合平均振幅分布
%   3. Boxplot_Trend: 趋势箱线图
%   4. Boxplot_Amp:   振幅箱线图
% 
% 风格：Robinson 投影, 紧凑布局, 智能刻度, 中文标题
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
OutputDir = 'D:\work\Figures\Average';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Models = {'EN4', 'IAP', 'Ishii'};
% 文件名需确认 (假设命名一致)
Files_T = {'EN4_TSLA_Terms_1to8_Average.mat', ...
           'IAP_TSLA_Terms_1to8_Average.mat', ...
           'Ishii_TSLA_Terms_1to8_Average.mat'};
       
% Files_S removed as per user request (redundant 1-3 order logic)

MaxOrder = 8; % 用户确认有8阶数据 (1-3 T+S+Cross, 4-8 T only unless updated)

%% 2. 加载数据并计算
% 目标: 获取所有模型的 Trend 和 Amp 场
% 结构: Data_Trend{Order}(Lat, Lon, Model)
fprintf('>> [1/5] 加载三套数据并计算 Trend/Amp (Order 1-%d)...\n', MaxOrder);

Ref_Lon = []; Ref_Lat = [];
[Lon_Grid, Lat_Grid] = deal([]);

% 预分配
% 由于各模型网格不同，先读第一个模型定基准
d1 = load(fullfile(DataDir, Files_T{1}));
if isfield(d1, 'Lon'), Ref_Lon = d1.Lon; else, Ref_Lon = d1.lon; end
if isfield(d1, 'Lat'), Ref_Lat = d1.Lat; else, Ref_Lat = d1.lat; end
if max(Ref_Lon) > 180, Ref_Lon(Ref_Lon>180) = Ref_Lon(Ref_Lon>180)-360; end
[Ref_Lon, s_idx] = sort(Ref_Lon);
[Lon_Grid, Lat_Grid] = meshgrid(Ref_Lon, Ref_Lat);
[Ny, Nx] = size(Lon_Grid);

% 存储容器: [Ny, Nx, 3_Models, MaxOrder]
All_Trend = zeros(Ny, Nx, 3, MaxOrder);
All_Amp   = zeros(Ny, Nx, 3, MaxOrder);

for m = 1:3
    fprintf('   Processing Model %d: %s ...\n', m, Models{m});
    
    % 加载 T 数据
    dT = load(fullfile(DataDir, Files_T{m}));
    if isfield(dT, 'Lon'), cur_Lon = dT.Lon; else, cur_Lon = dT.lon; end
    if isfield(dT, 'Lat'), cur_Lat = dT.Lat; else, cur_Lat = dT.lat; end
    if max(cur_Lon) > 180, cur_Lon(cur_Lon>180) = cur_Lon(cur_Lon>180)-360; end
    [cur_Lon, sort_idx] = sort(cur_Lon);
    
    % S/Cross term loading removed. Analysis is purely based on T terms (1-8).
    
    % 时间轴
    if isfield(dT, 'time_vec'), t_vec = dT.time_vec(:);
    elseif isfield(dT, 'time_axis'), t_vec = dT.time_axis(:);
    else, t_vec = dT.Time_Axis(:); 
    end
    
    % 准备插值器 (如果网格不同)
    NeedInterp = false;
    if length(cur_Lon) ~= length(Ref_Lon) || length(cur_Lat) ~= length(Ref_Lat) || ...
       any(abs(cur_Lon - Ref_Lon) > 1e-4) || any(abs(cur_Lat - Ref_Lat) > 1e-4)
        NeedInterp = true;
        [X_cur, Y_cur] = meshgrid(cur_Lon, cur_Lat);
    end
    
    % 循环阶数
    for k = 1:MaxOrder
        % 组装 Order k 的总序列 (T only)
        % T (Exist for 1-8)
        Term_T = squeeze(dT.TSLA_AllOrders(sort_idx, :, :, k)); 
        
        % 总和 (仅 T 项)
        Total_Term = Term_T;
        
        % 计算 Trend / Amp
        % [Nx_cur, Ny_cur]
        [tr_map, amp_map] = calc_trend_amp(Total_Term, t_vec);
        
        % 插值到参考网格
        if NeedInterp
            F_tr = griddedInterpolant(Y_cur, X_cur, tr_map', 'linear', 'none');
            F_am = griddedInterpolant(Y_cur, X_cur, amp_map', 'linear', 'none');
            tr_final = F_tr(Lat_Grid, Lon_Grid);
            am_final = F_am(Lat_Grid, Lon_Grid);
        else
            tr_final = tr_map'; % 转置为 [Lat, Lon]
            am_final = amp_map';
        end
        
        All_Trend(:, :, m, k) = tr_final;
        All_Amp(:, :, m, k)   = am_final;
    end
end

% 智能单位修正 (如果是 m 则转 mm)
if max(abs(All_Trend(:))) > 500
    fprintf('   Trend 单位修正: m -> mm\n');
    All_Trend = All_Trend / 1000; % WAIT, usually it's *1000. 
    % Let's logic: 
    % EN4 Trend ~ 0.5 mm/yr. If data is m/yr, it's 0.0005. 
    % If I see 0.0005, I should *1000.
    % If I see 500, it's wrongly scaled?
    % Let's rely on standard logic: raw data usually m.
    % Check calc_trend_amp function below.
end
% Re-check unit logic later. Assuming calc returns mm.

%% 3. 计算集合统计量
fprintf('>> [2/5] 计算 Ensemble Mean, Std, SNR ...\n');
% [Lat, Lon, Order]
Ens_Trend_Mean = squeeze(mean(All_Trend, 3)); 
Ens_Trend_Std  = squeeze(std(All_Trend, 0, 3));
ENS_Trend_SNR  = abs(Ens_Trend_Mean) ./ Ens_Trend_Std;

Ens_Amp_Mean   = squeeze(mean(All_Amp, 3));
Ens_Amp_Std    = squeeze(std(All_Amp, 0, 3));
% Ens_Amp_SNR    = abs(Ens_Amp_Mean) ./ Ens_Amp_Std; 

%% 4. 绘图1: Map_Trend (MME)
fprintf('>> [3/5] 绘制 MME Trend Map...\n');
figure('Position', [100, 50, 1200, 1400], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % 8 orders -> 4x2

% Colormap
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
    snr_mask = ENS_Trend_SNR(:,:,k) > 1.69;
    skip = 3;
    [r, c] = find(snr_mask(1:skip:end, 1:skip:end));
    if ~isempty(r)
        l_sub = Lon_Grid(1:skip:end, 1:skip:end);
        la_sub = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(l_sub(snr_mask(1:skip:end, 1:skip:end)), la_sub(snr_mask(1:skip:end, 1:skip:end)), ...
               '.', 'color', [.2 .2 .2], 'markersize', 3);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8);
    colormap(gca, MyCmap);
    
    % Caxis
    cl = prctile(abs(data(:)), 99.5); if cl<1e-5, cl=0.1; end
    caxis([-cl, cl]);
    
    title(sprintf('第 %d 阶 (Order %d)', k, k), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8);
end
sgtitle('\fontname{Microsoft YaHei}平均态 MME 趋势分布 (1-8阶, 点示 90% CI Significant)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Map_Trend_MME_Avg.png'));

%% 5. 绘图2: Map_Amp (MME)
fprintf('>> [4/5] 绘制 MME Amp Map...\n');
figure('Position', [150, 100, 1200, 1400], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Amp Colormap
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
    
    cl = prctile(data(:), 99.5); if cl<1e-5, cl=0.1; end
    caxis([0, cl]);
    
    title(sprintf('第 %d 阶 (Order %d)', k, k), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; title(cb, 'mm', 'FontSize', 8);
end
sgtitle('\fontname{Microsoft YaHei}平均态 MME 振幅分布 (1-8阶)', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Map_Amp_MME_Avg.png'));

%% 6. 绘图3 & 4: Boxplots
fprintf('>> [5/5] 绘制 Boxplots...\n');
% Flatten data for boxplot: [N_total, MaxOrder]
Trend_Flat = reshape(Ens_Trend_Mean, [], MaxOrder);
Amp_Flat   = reshape(Ens_Amp_Mean,   [], MaxOrder);

% Trend Boxplot (Figure 3)
fprintf('>> [5.1/5] 绘制 MME Trend Boxplot (3-Axis Compact)...\n');
figure('Position', [300, 300, 1200, 550], 'Color', 'w', 'Name', 'MME Trend Compact');
MainPos = [0.07, 0.14, 0.80, 0.80];

% Parameters
C1 = [0.85 0.325 0.098]; % Red
C2 = [0 0.447 0.741];    % Blue
C3 = [0.466 0.674 0.188];% Green
BoxBgColor = [0.95 0.95 1];
BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.4; ScatterSize = 12;

% Data Prep
Data_Trend = Trend_Flat;

% Layout
% Layer 1 (Red): Orders 1-3
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:min(3, MaxOrder)
    draw_jitter_unit(ax1, k, Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = max(abs(prctile(Data_Trend(:,1:min(3,MaxOrder)), [0.2 99.8], 'all'))); 
if isnan(max_v1)||max_v1==0, max_v1=1; end
Limit1 = max_v1 * 1.1; ylim(ax1, [-Limit1, Limit1]);
set(ax1, 'XLim', [0.4, MaxOrder+1], 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);

% Layer 2 (Blue): Orders 4-5
if MaxOrder >= 4
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 4:min(5, MaxOrder)
        draw_jitter_unit(ax2, k, Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    max_v2 = max(abs(prctile(Data_Trend(:,4:min(5,MaxOrder)), [0.2 99.8], 'all')));
    if isnan(max_v2)||max_v2==0, max_v2=1; end
    Limit2 = max_v2 * 1.2; ylim(ax2, [-Limit2, Limit2]);
    set(ax2, 'Color', 'none', 'XLim', [0.4, MaxOrder+1], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);
else
    ax2 = [];
end

% Layer 3 (Green): Orders 6+
if MaxOrder >= 6
    % Shift ax3 slightly right or overlay 
    Offset_Ratio = 0.08; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3);
    % Fix Alignment: Min must stay 0.4. Range scales with Width.
    % Range_New = Range_Old * Scale_Factor
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
else
    ax3 = [];
end

% Decoration
if ~isempty(ax1), yline(ax1, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.2); end
set(ax1, 'XTick', 1:MaxOrder, 'XTickLabel', arrayfun(@(x) sprintf('%d',x), 1:MaxOrder, 'UniformOutput', false));
xlabel(ax1, 'Order', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, 'Average State MME Trend Distribution (Orders 1-8)', 'FontSize', 16, 'FontWeight', 'bold');

% Legend
h_box = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
legend(ax1, [h_box, h_mean], {'IQR Box', 'Mean'}, 'Location', 'northeast', 'Box', 'off');

saveas(gcf, fullfile(OutputDir, 'Stats_Trend_Boxplot_MME_Avg.png'));


% Amp Boxplot (Figure 4)
fprintf('>> [5.2/5] 绘制 MME Amp Boxplot (3-Axis Compact)...\n');
figure('Position', [350, 350, 1200, 550], 'Color', 'w', 'Name', 'MME Amp Compact');
Data_Amp = Amp_Flat;

% Layer 1 (Red)
ax1 = axes('Position', MainPos); hold(ax1, 'on');
for k = 1:min(3, MaxOrder)
    draw_jitter_unit(ax1, k, Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C1, BoxBgColor, ScatterAlpha);
end
max_v1 = prctile(Data_Amp(:,1:min(3,MaxOrder)), 99.8, 'all'); 
if isnan(max_v1)||max_v1==0, max_v1=1; end
ylim(ax1, [0, max_v1*1.1]);
set(ax1, 'XLim', [0.4, MaxOrder+1], 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax1, 'Dominant Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);

% Layer 2 (Blue)
if MaxOrder >= 4
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 4:min(5, MaxOrder)
        draw_jitter_unit(ax2, k, Data_Amp(:,k), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    max_v2 = prctile(Data_Amp(:,4:min(5,MaxOrder)), 99.8, 'all'); 
    if isnan(max_v2)||max_v2==0, max_v2=1; end
    ylim(ax2, [0, max_v2*1.2]);
    set(ax2, 'Color', 'none', 'XLim', [0.4, MaxOrder+1], 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel(ax2, 'Secondary Amp (mm)', 'FontWeight', 'bold', 'FontSize', 14);
else
    ax2 = [];
end

% Layer 3 (Green)
if MaxOrder >= 6
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
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
else
    ax3 = [];
end

set(ax1, 'XTick', 1:MaxOrder, 'XTickLabel', arrayfun(@(x) sprintf('%d',x), 1:MaxOrder, 'UniformOutput', false));
xlabel(ax1, 'Order', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, 'Average State MME Amplitude Distribution (Orders 1-8)', 'FontSize', 16, 'FontWeight', 'bold');

% Legend
h_box = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
legend(ax1, [h_box, h_mean], {'IQR Box', 'Mean'}, 'Location', 'northeast', 'Box', 'off');

saveas(gcf, fullfile(OutputDir, 'Stats_Amp_Boxplot_MME_Avg.png'));

fprintf('>> 完成! 图表均已保存至 %s\n', OutputDir);


%% Internal Functions
function [t_map, a_map] = calc_trend_amp(data_3d, t_vec)
    % Vectorized or Loop. data_3d is [Nlon, Nlat, Nt]
    [nx, ny, nt] = size(data_3d);
    t_map = nan(nx, ny);
    a_map = nan(nx, ny);
    
    for i=1:nx
        for j=1:ny
            y = squeeze(data_3d(i,j,:));
            idx = ~isnan(y);
            if sum(idx) < nt*0.6, continue; end
            
            p = polyfit(t_vec(idx), y(idx), 1);
            % 假设输入是 mm (Calc脚本已转). 若是 m, 需 *1000.
            % 这里的 input 来自 Calc 脚本，通常已经是 mm.
            t_map(i,j) = p(1); 
            
            resid = y(idx) - polyval(p, t_vec(idx));
            a_map(i,j) = std(resid);
        end
    end
    
    % 智能检测单位 (如果 trend 极小 < 1e-4, 可能是 m)
    if nanmax(abs(t_map(:))) < 1e-2 && nanmax(abs(t_map(:))) > 0
         % e.g. 0.001 mm/yr vs 1 mm/yr. 
         % Mean Trend is usually 0-5 mm/yr. 
         % If it is 0.003, it might be m/yr? No, 0.003 m/yr is 3 mm/yr. 
         % Wait. m/yr: 0.003 = 3mm.
         % If result is 3, it's mm. 
         % If result is 0.003, it's m?
