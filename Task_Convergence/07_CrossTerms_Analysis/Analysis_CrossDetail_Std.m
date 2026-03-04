%% Analysis_CrossDetail_Std.m
% =========================================================================
% 功能：标准态混合项 (2-8阶, 28项) 详细分析
% 内容：
%   1. Trend & Amp 空间分布 (Total 28 terms, Split into 3 Pages PDF)
%   2. Boxplot 对比
%   3. Decay Curve (衰减曲线)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\CrossDetail_Std';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

File = fullfile(DataDir, 'EN4_CrossDetail_Std.mat');
if ~exist(File, 'file'), error('请先运行 Calc_CrossDetail_Std.m'); end

fprintf('>> [1/5] Loading Data...\n');
D = load(File);

% 整理变量列表 (按阶排序)
VarList = {}; Titles = {}; Orders = [];
for n = 2:8
    for k = 1:n-1
        VarList{end+1} = sprintf('Cross_T%dS%d', n-k, k);
        Titles{end+1} = sprintf('%d阶: T^%dS^%d', n, n-k, k);
        Orders(end+1) = n;
    end
end
NumVars = length(VarList); % 28

% 网格
Lon = D.Lon; Lat = D.Lat; t_vec = D.time_vec(:);
if max(Lon)>180, Lon(Lon>180)=Lon(Lon>180)-360; end
[Lon, sort_idx] = sort(Lon);
[LonGrid, LatGrid] = meshgrid(Lon, Lat);

%% 2. 计算统计量
fprintf('>> [2/5] 正在计算统计量...\n');
Trend_Stats = struct('Max', [], 'Mean', [], 'P99', [], 'Order', []);
Amp_Stats = struct('Max', [], 'Mean', [], 'P99', [], 'Order', []);

Trend_Maps = cell(1, NumVars);
Amp_Maps = cell(1, NumVars);

calc_fun = @(d) compute_stats_vectorized(d, t_vec);

for k = 1:NumVars
    name = VarList{k};
    if ~isfield(D, name)
        Trend_Maps{k} = nan(size(LatGrid')); Amp_Maps{k} = nan(size(LatGrid'));
        continue; 
    end
    data = D.(name)(sort_idx, :, :);
    [Tr, Am, ~] = calc_fun(data);
    Trend_Maps{k} = Tr;
    Amp_Maps{k} = Am;
    
    % Stats for Decay Curve
    v_tr = abs(Tr(:)); v_tr(isnan(v_tr))=[];
    v_am = Am(:); v_am(isnan(v_am))=[];
    
    if isempty(v_tr), v_tr = 0; end
    if isempty(v_am), v_am = 0; end
    
    Trend_Stats(k).Max = max(v_tr);
    Trend_Stats(k).Mean = mean(v_tr);
    Trend_Stats(k).P99 = prctile(v_tr, 99);
    Trend_Stats(k).Order = Orders(k);
    
    Amp_Stats(k).Max = max(v_am);
    Amp_Stats(k).Mean = mean(v_am);
    Amp_Stats(k).P99 = prctile(v_am, 99);
    Amp_Stats(k).Order = Orders(k);
end

%% 3. 绘图 - Trend Maps (分页 PDF)
fprintf('>> [3/5] 绘制 Trend 分布 (分页 PDF)...\n');

% 分组定义
Groups = {2:5, 6:7, 8};
GroupNames = {'2-5阶', '6-7阶', '8阶'};
GroupLayouts = {[3, 4], [3, 4], [2, 4]}; % Rows, Cols
PDF_Name = fullfile(OutputDir, 'Trend_Maps_Std_All.pdf');
if exist(PDF_Name, 'file'), delete(PDF_Name); end

% Blue-White-Red Cmap
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx_c = linspace(0, 1, 256);
Cmap = [interp1(ColorLocs, ColorPoints(:,1), xx_c, 'pchip')', ...
        interp1(ColorLocs, ColorPoints(:,2), xx_c, 'pchip')', ...
        interp1(ColorLocs, ColorPoints(:,3), xx_c, 'pchip')'];

for g = 1:length(Groups)
    CurrentOrders = Groups{g};
    Mask = ismember(Orders, CurrentOrders);
    Indices = find(Mask);
    
    if isempty(Indices), continue; end
    
    figure('Position', [50, 50, 1600, 900], 'Color', 'w');
    layout = GroupLayouts{g};
    t = tiledlayout(layout(1), layout(2), 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for k = Indices
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        m_pcolor(LonGrid, LatGrid, Trend_Maps{k}); shading flat;
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
            'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
        colormap(gca, Cmap);
        title(Titles{k}, 'FontSize', 10, 'Interpreter', 'tex', 'FontWeight', 'bold');
        
        v = Trend_Maps{k}(:); clim = prctile(abs(v), 99);
        if isnan(clim) || clim < 1e-4, clim = 1e-4; end
        caxis([-clim, clim]);
        cb = colorbar; title(cb, 'mm/yr', 'FontSize', 7);
    end
    
    sgtitle(sprintf('标准态混合项趋势 (%s)', GroupNames{g}), 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(gcf, PDF_Name, 'Append', true);
    % Keep open
end


%% 4. 绘图 - Amp Maps (分页 PDF)
fprintf('>> [4/5] 绘制 Amp 分布 (分页 PDF)...\n');

PDF_Name_Amp = fullfile(OutputDir, 'Amp_Maps_Std_All.pdf');
if exist(PDF_Name_Amp, 'file'), delete(PDF_Name_Amp); end

% Red Seq Cmap
MyRedSeq = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx = [0 0.1 0.4 0.8 1]; yy = linspace(0,1,256);
CmapSeq = [interp1(xx,MyRedSeq(:,1),yy)', interp1(xx,MyRedSeq(:,2),yy)', interp1(xx,MyRedSeq(:,3),yy)'];

for g = 1:length(Groups)
    CurrentOrders = Groups{g};
    Mask = ismember(Orders, CurrentOrders);
    Indices = find(Mask);
    
    if isempty(Indices), continue; end
    
    figure('Position', [50, 50, 1600, 900], 'Color', 'w');
    layout = GroupLayouts{g};
    t = tiledlayout(layout(1), layout(2), 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for k = Indices
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        m_pcolor(LonGrid, LatGrid, Amp_Maps{k}); shading flat;
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
            'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
        colormap(gca, CmapSeq);
        title(Titles{k}, 'FontSize', 10, 'Interpreter', 'tex', 'FontWeight', 'bold');
        
        v = Amp_Maps{k}(:); clim = prctile(v, 99);
        if isnan(clim) || clim < 1e-4, clim = 1e-4; end
        caxis([0, clim]);
        cb = colorbar; title(cb, 'mm', 'FontSize', 7);
    end
    
    sgtitle(sprintf('标准态混合项振幅 (%s)', GroupNames{g}), 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(gcf, PDF_Name_Amp, 'Append', true);
    % Keep open
end

%% 5. 衰减曲线 (Decay Curve - 核心图表)
fprintf('>> [5/5] 绘制衰减曲线...\n');

OrderList = 2:8;
MaxTrendPerOrder = zeros(1, 7);
MeanTrendPerOrder = zeros(1, 7);
MaxAmpPerOrder = zeros(1, 7);
P99TrendPerOrder = zeros(1, 7);

for i = 1:length(OrderList)
    ord = OrderList(i);
    idx = find([Trend_Stats.Order] == ord);
    
    MaxTrendPerOrder(i) = max([Trend_Stats(idx).Max]);
    MeanTrendPerOrder(i) = mean([Trend_Stats(idx).Mean]);
    P99TrendPerOrder(i) = max([Trend_Stats(idx).P99]); % Use Max of P99 to be safe
    MaxAmpPerOrder(i) = max([Amp_Stats(idx).Max]);
end

figure('Position', [100, 100, 800, 400], 'Color', 'w');
subplot(1,2,1);
semilogy(OrderList, MaxTrendPerOrder, 'r-o', 'LineWidth', 2, 'DisplayName', 'Max Trend');
hold on;
semilogy(OrderList, P99TrendPerOrder, 'b-s', 'LineWidth', 1.5, 'DisplayName', 'P99 Trend');
semilogy(OrderList, MeanTrendPerOrder, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'Mean Trend');
grid on; xlabel('Order'); ylabel('Trend (mm/yr)'); legend('Location', 'ne');
    title('标准态混合项趋势随阶数衰减 (Log Scale)');
    
    subplot(1,2,2);
    semilogy(OrderList, MaxAmpPerOrder, 'm-d', 'LineWidth', 2);
    grid on; xlabel('Order'); ylabel('Amplitude (mm)');
    title('标准态混合项最大振幅随阶数衰减');
    
    sgtitle('标准态高阶混合项收敛性检验');
saveas(gcf, fullfile(OutputDir, 'Decay_Curve_Std.png'));

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
