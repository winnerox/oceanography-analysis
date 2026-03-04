%% Analysis_CrossDetail_Avg.m
% =========================================================================
% 功能：平均态混合项 (2-5阶, 10项) 详细分析
% 内容：
%   1. Trend & Amp 空间分布 (10 terms)
%   2. Boxplot 对比
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\CrossDetail_Avg';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

File = fullfile(DataDir, 'EN4_CrossDetail_Avg.mat');
if ~exist(File, 'file'), error('请先运行 Calc_CrossDetail_Avg.m'); end

fprintf('>> [1/4] Loading Data...\n');
D = load(File);

% 整理变量列表 (按阶排序)
VarList = {};
Titles = {};
% 2阶: T1S1
VarList{end+1} = 'Cross_T1S1'; Titles{end+1} = '2阶: T^{1}S^{1}';
% 3阶: T2S1, T1S2
VarList{end+1} = 'Cross_T2S1'; Titles{end+1} = '3阶: T^{2}S^{1}';
VarList{end+1} = 'Cross_T1S2'; Titles{end+1} = '3阶: T^{1}S^{2}';
% 4阶: T3S1, T2S2, T1S3
VarList{end+1} = 'Cross_T3S1'; Titles{end+1} = '4阶: T^{3}S^{1}';
VarList{end+1} = 'Cross_T2S2'; Titles{end+1} = '4阶: T^{2}S^{2}';
VarList{end+1} = 'Cross_T1S3'; Titles{end+1} = '4阶: T^{1}S^{3}';
% 5阶: T4S1..T1S4
VarList{end+1} = 'Cross_T4S1'; Titles{end+1} = '5阶: T^{4}S^{1}';
VarList{end+1} = 'Cross_T3S2'; Titles{end+1} = '5阶: T^{3}S^{2}';
VarList{end+1} = 'Cross_T2S3'; Titles{end+1} = '5阶: T^{2}S^{3}';
VarList{end+1} = 'Cross_T1S4'; Titles{end+1} = '5阶: T^{1}S^{4}';

NumVars = length(VarList); % 10

% Grid
Lon = D.Lon; Lat = D.Lat; t_vec = D.time_vec(:);
if max(Lon)>180, Lon(Lon>180)=Lon(Lon>180)-360; end
[Lon, sort_idx] = sort(Lon);
[LonGrid, LatGrid] = meshgrid(Lon, Lat);

%% 2. 计算统计量
fprintf('>> [2/4]正在计算统计量...\n');
Trend_Maps = cell(1, NumVars);
Amp_Maps = cell(1, NumVars);

calc_fun = @(d) compute_stats_vectorized(d, t_vec);

for k = 1:NumVars
    name = VarList{k};
    if ~isfield(D, name), continue; end
    data = D.(name)(sort_idx, :, :);
    [Tr, Am, ~] = calc_fun(data);
    Trend_Maps{k} = Tr;
    Amp_Maps{k} = Am;
end

%% 3. 绘制分布图
fprintf('>> [3/4] 正在绘图...\n');

% Trend Map (2x5)
figure('Position', [100, 100, 1800, 600], 'Color', 'w');
t = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

% Blue-White-Red Cmap (Same as Map_Trend_EN4.m)
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
Cmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
        interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
        interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];

for k = 1:NumVars
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(LonGrid, LatGrid, Trend_Maps{k}); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
    colormap(gca, Cmap);
    title(Titles{k}, 'FontSize', 10, 'FontWeight', 'bold');
    
    % Caxis: Auto 98%
    v = Trend_Maps{k}(:); clim = prctile(abs(v), 98);
    if isnan(clim) || clim < 1e-3, clim = 1e-3; end
    caxis([-clim, clim]);
    cb = colorbar; title(cb, 'mm/yr', 'FontSize', 8);
end
sgtitle('平均态混合项趋势分布 (2-5阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Trend_Maps_Avg_Detail.png'));

% Amp Map
figure('Position', [100, 100, 1800, 600], 'Color', 'w');
t = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
% Red Seq Cmap (Same as Map_Amp_EN4.m)
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx = [0 0.1 0.4 0.8 1]; yy = linspace(0,1,256);
CmapSeq = [interp1(xx,MyRedSeq(:,1),yy)', interp1(xx,MyRedSeq(:,2),yy)', interp1(xx,MyRedSeq(:,3),yy)'];

for k = 1:NumVars
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(LonGrid, LatGrid, Amp_Maps{k}); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
    colormap(gca, CmapSeq);
    title(Titles{k}, 'FontSize', 10, 'FontWeight', 'bold');
    
    v = Amp_Maps{k}(:); clim = prctile(v, 98);
    if isnan(clim) || clim < 1e-3, clim = 1e-3; end
    caxis([0, clim]);
    cb = colorbar; title(cb, 'mm', 'FontSize', 8);
end
sgtitle('平均态混合项振幅分布 (2-5阶)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(OutputDir, 'Amp_Maps_Avg_Detail.png'));

%% 4. 绘制箱线图
fprintf('>> [4/4] 绘制箱线图...\n');
figure('Position', [100, 100, 1200, 500], 'Color', 'w');
subplot(1,2,1); hold on;
for k = 1:NumVars
    data = Trend_Maps{k}(:); data(isnan(data))=[];
    if isempty(data), continue; end
    boxplot(data, 'Positions', k, 'Widths', 0.6, 'Symbol', '');
end
xlim([0, NumVars+1]); xticks(1:NumVars); xticklabels(Titles); xtickangle(45);
ylabel('Trend (mm/yr)'); title('平均态混合项趋势分布详情'); grid on;

subplot(1,2,2); hold on;
for k = 1:NumVars
    data = Amp_Maps{k}(:); data(isnan(data))=[];
    if isempty(data), continue; end
    boxplot(data, 'Positions', k, 'Widths', 0.6, 'Symbol', '');
end
xlim([0, NumVars+1]); xticks(1:NumVars); xticklabels(Titles); xtickangle(45);
ylabel('Amp (mm)'); title('平均态混合项振幅分布详情'); grid on;

sgtitle('平均态混合项统计箱线图');
saveas(gcf, fullfile(OutputDir, 'Boxplot_Avg_Detail.png'));

fprintf('>> Done!\n');

function [trend, amp, sig] = compute_stats_vectorized(data_3d, t_vec)
    [Nx, Ny, Nt] = size(data_3d);
    t_vec = t_vec(:); 
    trend = nan(Nx, Ny); amp = nan(Nx, Ny); sig = nan(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            y = squeeze(data_3d(i, j, :));
            mask = ~isnan(y);
            if sum(mask) < Nt * 0.7, continue; end
            x = t_vec(mask); y = y(mask);
            p = polyfit(x, y, 1);
            trend(i, j) = p(1);
            amp(i, j) = std(y - polyval(p, x));
        end
    end
    trend = trend'; amp = amp'; sig = sig';
end
