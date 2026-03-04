function Plot_Hybrid_Trend_Combined(dataset_name, state)
% 功能: 绘制混合项趋势空间分布图（无箱线图版），自适应多阶数量
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 绘制混合项趋势空间分布图 ==========\n');

%% 1. 加载趋势结果文件
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
trend_file = fullfile(trend_dir, sprintf('%s_%s_Hybrid_Trends.mat', dataset_name, state));

if ~exist(trend_file, 'file')
    error('❌ 找不到混合项趋势结果文件: %s', trend_file);
end

fprintf('>> 加载混合项趋势数据: %s ...\n', trend_file);
load(trend_file);

%% 2. 准备数据
if ~exist('trend_Hybrid', 'var')
    error('❌ 数据文件中未找到 trend_Hybrid 变量，请检查 Calculate_Hybrid_Trends 的输出。');
end
trend_data = trend_Hybrid;
sig_data = sig_Hybrid;

% 经纬度处理
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end

% 确保Lon和Lat是向量
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180] 完美对齐
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    for i = 1:size(trend_data, 3)
        trend_data(:,:,i) = trend_data(:, sort_idx, i);
        sig_data(:,:,i) = sig_data(:, sort_idx, i);
    end
end
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

% 转置使得画图时不报错
if size(trend_data, 1) ~= length(Lat)
    for i = 1:size(trend_data, 3)
        temp_t(:,:,i) = trend_data(:,:,i)';
        temp_s(:,:,i) = sig_data(:,:,i)';
    end
    trend_data = temp_t;
    sig_data = temp_s;
end

%% 3. 绘制空间分布图 (顶刊 RdBu 色棒)
norder = size(trend_data, 3);
fprintf('>> 检测到 %d 个混合项阶数，准备动态分配画布布局...\n', norder);

% 自适应决定子图布局 (rows, cols)
if norder <= 3
    rows = 1; cols = norder;
elseif norder == 4
    rows = 2; cols = 2;
elseif norder <= 6
    rows = 2; cols = 3;
elseif norder <= 8
    rows = 4; cols = 2;
else
    cols = ceil(sqrt(norder));
    rows = ceil(norder / cols);
end

% 动态调整画板大小
fig_width = min(400 * cols, 1600);
fig_height = min(350 * rows + 100, 1200);

figure('Position', [100, 100, fig_width, fig_height], 'Color', 'w', 'Name', sprintf('Hybrid Trend Spatial Maps - %s %s', dataset_name, state));
t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘 (深蓝-浅蓝-纯白-浅红-深红)
RdBu_ColorPoints = [
    0.1 0.2 0.5;
    0.6 0.8 0.9;
    1.0 1.0 1.0;
    1.0 0.7 0.7;
    0.6 0.1 0.1
];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; 
xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:norder
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = trend_data(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % 显著性打点
    sig_mask = sig_data(:,:,term_idx);
    skip = 3; 
    [row, col] = find(sig_mask(1:skip:end, 1:skip:end) == 1);
    if ~isempty(row)
        sig_lon = Lon_Grid(1:skip:end, 1:skip:end);
        sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==1), ...
               sig_lat(sig_mask(1:skip:end, 1:skip:end)==1), ...
               'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, MyDivergentCmap);
    
    % 计算适当的色标上限
    raw_max = prctile(abs(data(:)), 99.8);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([-final_limit, final_limit]); 
    
    title(sprintf('Hybrid Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) Hybrid Terms Trend Distribution (p < 0.05)', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 绘图完成！\n');
end
