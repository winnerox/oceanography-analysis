%% 振幅空间分布图
%% Step40_Amp_4x2_Compact_HighContrast.m
% 功能: 振幅空间分布图 (完美复刻趋势图样式：紧凑布局 + 强显色)
% 改进:
%   1. 【单位置顶】将单位 'mm' 移至 Colorbar 顶部，与趋势图保持一致。
%   2. 【显色增强】分位数阈值降至 95% (原98%)，解决"太白"问题，让弱信号显色。
%   3. 【智能刻度】保持整洁的整数或 .5 结尾刻度。
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

% 计算 Amp
Map_Amp = nan(nlat, nlon, 8);
t_vec = time_axis(:);
fprintf('>> 计算振幅分布 (mm)...\n');
for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    temp = nan(nlon, nlat);
    for i=1:nlon, for j=1:nlat
        ts = squeeze(slice(i,j,:));
        idx_valid = ~isnan(ts);
        if sum(idx_valid) < length(ts) * 0.5, continue; end
        
        p = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
        detrended_ts = ts(idx_valid) - polyval(p, t_vec(idx_valid));
        temp(i,j) = std(detrended_ts) * 1000;
    end; end
    Map_Amp(:,:,n) = temp';
end

%% ==========================================
%% 绘图 (4x2 + 紧凑单位 + 强显色)
%% ==========================================
figure('Position', [100, 50, 1200, 1400], 'Color', 'w'); 
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% 振幅专用色标 (白->红)
MyRed = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx = [0 0.1 0.4 0.8 1]; yy = linspace(0,1,256);
MySeqCmap = [interp1(xx,MyRed(:,1),yy)', interp1(xx,MyRed(:,2),yy)', interp1(xx,MyRed(:,3),yy)'];

for term_idx = 1:8
    nexttile;
    
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = Map_Amp(:,:,term_idx);
    
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
    
    colormap(gca, MySeqCmap);
    
    % === 智能取整算法 (振幅版) ===
    % 1. 使用 95% 分位数 (与趋势图保持一致的强显色策略)
    raw_max = prctile(data(:), 95.0); 
    if raw_max < 1e-12, raw_max = 1e-12; end
    
    % 2. 计算数量级
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    
    % 3. 向上取整 (增加细分档位，避免颜色过浅)
    if fraction <= 1
        nice_base = 1;
    elseif fraction <= 1.25
        nice_base = 1.25;
    elseif fraction <= 1.5
        nice_base = 1.5;
    elseif fraction <= 2
        nice_base = 2;
    elseif fraction <= 2.5
        nice_base = 2.5;
    elseif fraction <= 3
        nice_base = 3;
    elseif fraction <= 4
        nice_base = 4;
    elseif fraction <= 5
        nice_base = 5;
    elseif fraction <= 6
        nice_base = 6;
    elseif fraction <= 8
        nice_base = 8;
    else
        nice_base = 10;
    end
    
    final_limit = nice_base * 10^exponent;
    
    % 设置范围 [0, Max]
    caxis([0, final_limit]); 
    
    % 标题
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % === 【核心修改】紧凑单位布局 ===
    cb = colorbar; 
    cb.FontSize = 8;
    
    % 将单位放在顶部
    title(cb, 'mm', 'FontSize', 8, 'FontName', 'Times New Roman');
    
    % 强制 5 个刻度: [0, 0.25L, 0.5L, 0.75L, L]
    cb.Ticks = linspace(0, final_limit, 5);
    
    try cb.Ruler.Exponent = 0; catch; end
end

sgtitle('\fontname{Microsoft YaHei}1-8 阶独立项振幅空间分布', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(gcf, 'PaperPositionMode', 'auto'); 

fprintf('>> 振幅图绘图完成！已完全对齐趋势图（紧凑布局+强显色）。\n');