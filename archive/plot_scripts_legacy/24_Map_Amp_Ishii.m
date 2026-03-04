%% Step21_Ishii_Map_Amp.m
% =========================================================================
% 功能: Ishii 数据振幅空间分布图 (8阶, 4x2布局)
% =========================================================================
clear; clc; close all;

%% 0. 配置
OutputDir = 'D:\work\Figures';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 1. 加载数据
DataDir = 'D:\work\MAT_Data';
files = {'Ishii_TSLA_Terms_1to8_Average.mat', 'Ishii_TSLA_Result.mat'};
input_file = '';
for i = 1:length(files)
    fp = fullfile(DataDir, files{i});
    if exist(fp, 'file'), input_file = fp; break; end
end
if isempty(input_file), error('❌ 找不到 Ishii 数据文件'); end
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
            temp(i,j) = std(detrended_ts) * 1000;
        end
    end
    Map_Amp(:,:,n) = temp';
end

%% 绘图
fprintf('>> 绘制 Ishii 振幅空间分布图...\n');
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

sgtitle('\fontname{Microsoft YaHei}Ishii 1-8 阶振幅空间分布', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

% 保存
OutFile = fullfile(OutputDir, 'Ishii_Amp_Map_8Orders.jpg');
exportgraphics(gcf, OutFile, 'Resolution', 300);
fprintf('>> 图片已保存: %s\n', OutFile);
