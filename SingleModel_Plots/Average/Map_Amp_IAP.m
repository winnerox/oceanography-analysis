%% 振幅空间分布图 [Average Version]
%% Map_Amp_IAP.m
% 功能: IAP 数据振幅空间分布图 (8阶, 4x2布局)
clear; clc; close all;

%% 0. 配置
OutputDir = 'D:\work\Figures';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 1. 加载数据
files = {'IAP_TSLA_Terms_1to8_Average.mat'};
input_file = '';
DataDir = 'D:\work\MAT_Data';
for i = 1:length(files)
    if exist(files{i}, 'file')
        input_file = files{i};
    elseif exist(fullfile(DataDir, files{i}), 'file')
        input_file = fullfile(DataDir, files{i});
    end
    if ~isempty(input_file), break; end
end
if isempty(input_file), error('❌ 找不到 IAP 数据文件: %s', files{1}); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

if ~exist('time_axis', 'var') && exist('Time_Axis', 'var')
    time_axis = Time_Axis;
end

% 经纬度修正
if max(lon) > 180, lon(lon > 180) = lon(lon > 180) - 360; end
[lon, sort_idx] = sort(lon);
TSLA_AllOrders = TSLA_AllOrders(sort_idx, :, :, :);
[nlon, nlat, ntime, max_order] = size(TSLA_AllOrders);
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);
t_vec = time_axis(:);

% 计算振幅
fprintf('>> 计算振幅...\n');
Map_Amp = nan(nlat, nlon, 8);
for n = 1:min(8, max_order)
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    temp = nan(nlon, nlat);
    for i=1:nlon
        for j=1:nlat
            ts = squeeze(slice(i,j,:));
            idx_valid = ~isnan(ts);
            if sum(idx_valid) < length(ts) * 0.5, continue; end
            p = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
            detrended_ts = ts(idx_valid) - polyval(p, t_vec(idx_valid));
            
            % 统一处理：默认视为 m -> mm (*1000)
            % 后续在循环外进行全局智能修正
            temp(i,j) = std(detrended_ts) * 1000;
        end
    end
    Map_Amp(:,:,n) = temp';
end

% === 全局智能单位修正 ===
% 如果转换后数值过大 (例如 > 500 mm)，说明原数据本身就是 mm，不需要 *1000
if max(abs(Map_Amp(:))) > 500
    fprintf('>> 检测到数值过大 (Max=%.1f), 推测原数据为 mm, 自动回调 /1000\n', max(abs(Map_Amp(:))));
    Map_Amp = Map_Amp / 1000;
end

%% 绘图
fprintf('>> 绘制 IAP 振幅空间分布图...\n');
figure('Position', [100, 50, 1200, 1400], 'Color', 'w'); 
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% 配色 (白->红)
MyRed = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx = [0 0.1 0.4 0.8 1]; yy = linspace(0,1,256);
MySeqCmap = [interp1(xx,MyRed(:,1),yy)', interp1(xx,MyRed(:,2),yy)', interp1(xx,MyRed(:,3),yy)'];

for term_idx = 1:8
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = Map_Amp(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8, 'xticklabels', [], 'yticklabels', []);
    
    colormap(gca, MySeqCmap);
    
    raw_max = prctile(data(:), 95.0);
    if raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 1, nice_base = 1;
    elseif fraction <= 2, nice_base = 2;
    elseif fraction <= 5, nice_base = 5;
    else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([0, final_limit]); 
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    cb = colorbar; cb.FontSize = 8;
    title(cb, 'mm', 'FontSize', 8); 
    cb.Ticks = linspace(0, final_limit, 5);
    try cb.Ruler.Exponent = 0; catch; end
end

sgtitle('\fontname{Microsoft YaHei}IAP (Average) 1-8 阶振幅空间分布', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 图片已绘制完成。\n');
