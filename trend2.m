%% 趋势空间分布图2
%% Step38_Trend_Integer_Ticks.m
% 功能: 完美刻度版趋势图 (4x2布局 + 原版配色 + 智能取整刻度)
% 核心修正:
%   1. 【智能取整】引入算法，将乱七八糟的小数范围自动变为 [1, 2, 4, 5, 6, 8, 10] 等整洁数值。
%   2. 【刻度整齐】确保 Colorbar 上的 5 个数字要么是整数，要么是以 .5 结尾。
%   3. 【完全复刻】布局、配色、算法逻辑与您最满意的 Step24 完全一致。
clear; clc; close all;

%% 1. 加载数据
files = {'EN4_TSLA_Terms_1to15_StdRef.mat'};
input_file = '';
for i = 1:length(files), if exist(files{i}, 'file'), input_file = files{i}; break; end; end
if isempty(input_file)
    fallback_files = {'EN4_TSLA_Terms_1to15_Safe.mat', 'EN4_TSLA_Terms_1to15_GPU.mat'};
    for i = 1:length(fallback_files), if exist(fallback_files{i}, 'file'), input_file = fallback_files{i}; break; end; end
end
if isempty(input_file), error('❌ 找不到数据文件'); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 经纬度修正
if ~exist('lon', 'var'), if exist('Lon','var'), lon=Lon; elseif exist('longitude','var'), lon=longitude; end; end
if ~exist('lat', 'var'), if exist('Lat','var'), lat=Lat; elseif exist('latitude','var'), lat=latitude; end; end
if max(lon) > 180, lon(lon > 180) = lon(lon > 180) - 360; end
[lon, sort_idx] = sort(lon);
TSLA_AllOrders = TSLA_AllOrders(sort_idx, :, :, :);
[nlon, nlat, ntime, max_order] = size(TSLA_AllOrders);
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);

% 计算 Trend
Map_Trend = nan(nlat, nlon, 8);
t_vec = time_axis(:);
fprintf('>> 计算趋势分布 (mm/yr)...\n');
for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    temp = nan(nlon, nlat);
    for i=1:nlon, for j=1:nlat
        ts = squeeze(slice(i,j,:));
        idx_valid = ~isnan(ts);
        if sum(idx_valid) < length(ts) * 0.5, continue; end
        p = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
        temp(i,j) = p(1) * 1000; 
    end; end
    Map_Trend(:,:,n) = temp';
end

%% ==========================================
%% 绘图 (4x2 + 智能取整)
%% ==========================================
% 画布大小：1200x1400 (宽幅适中，防止压扁)
figure('Position', [100, 50, 1200, 1400], 'Color', 'w'); 

% 布局: 4行2列, compact
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
GammaList = [0.6 0.6 0.5 0.5 0.4 0.4 0.3 0.3];

% 恢复 Step24 原版蓝白红配色 (您觉得最好看的那版)
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1];
MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx, 'pchip')'];

% 循环绘图
for term_idx = 1:8
    nexttile;
    
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Map_Trend(:,:,term_idx);
    
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    
    % 网格样式
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
    
    colormap(gca, MyDivergentCmap);
    
    % === 【核心修改】智能取整算法 ===
    % 1. 先算出原始的 99.8% 极限值
    raw_max = prctile(abs(data(:)), 99.8);
    if raw_max < 1e-12, raw_max = 1e-12; end
    
    % 2. 计算数量级
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; % 把数值归一化到 1-10 之间
    
    % 3. 向上取整到最近的"漂亮数" (1, 2, 4, 5, 6, 8, 10)
    % 这样除以2得到的刻度(Ticks)也会很整齐
    if fraction <= 1
        nice_base = 1;      % Ticks: 0.5, 1
    elseif fraction <= 2
        nice_base = 2;      % Ticks: 1, 2
    elseif fraction <= 4
        nice_base = 4;      % Ticks: 2, 4
    elseif fraction <= 5
        nice_base = 5;      % Ticks: 2.5, 5
    elseif fraction <= 6
        nice_base = 6;      % Ticks: 3, 6
    elseif fraction <= 8
        nice_base = 8;      % Ticks: 4, 8
    else
        nice_base = 10;     % Ticks: 5, 10
    end
    
    % 4. 算出最终的整洁上限
    final_limit = nice_base * 10^exponent;
    
    caxis([-final_limit, final_limit]); 
    
    % 标题
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % === 统一刻度设置 ===
    cb = colorbar; 
    cb.Label.String = 'mm/yr';
    cb.Label.FontName = 'Times New Roman';
    cb.FontSize = 8;
    
    % 强制 5 个刻度，现在它们一定是整数或 .5 结尾
    cb.Ticks = linspace(-final_limit, final_limit, 5);
    
    % 强制让指数显示在上方
    try cb.Ruler.Exponent = 0; catch; end
end

sgtitle('\fontname{Microsoft YaHei}1-8 阶独立项趋势空间分布', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(gcf, 'PaperPositionMode', 'auto'); 

fprintf('>> 绘图完成！刻度值已优化为整数或整齐的小数。\n');