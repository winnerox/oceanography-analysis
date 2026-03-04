clc
clear
close all

load('total_info_Amp.mat')
load('total_info_Trend.mat')
load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\EN4_c13\EN4_c13_all_new.mat','lon','lat')


for k = 1 : length(total_info_Trend)
    total_info_Trend(k).Trend_2000a = cat(3,total_info_Trend(k).EN4_Trend_2000a, total_info_Trend(k).IAP_Trend_2000a, total_info_Trend(k).Ishii_Trend_2000a);
    total_info_Trend(k).Trend_2000b = cat(3,total_info_Trend(k).EN4_Trend_2000b, total_info_Trend(k).IAP_Trend_2000b, total_info_Trend(k).Ishii_Trend_2000b);
    total_info_Trend(k).Trend_all = cat(3,total_info_Trend(k).EN4_Trend_all, total_info_Trend(k).IAP_Trend_all, total_info_Trend(k).Ishii_Trend_all);
end

Trend_data = {total_info_Trend.Trend_2000a};

% ========== 第二部分：绘制所有子图 ==========
% 设置颜色映射
mycolor = othercolor('RdYlBu10',16);
mycolor = flipud(mycolor);

% 创建图形窗口（调整大小以容纳紧凑的子图）
figure('Position', [100, 100, 1400, 900]);

% 获取经纬度范围（假设所有层使用相同的经纬度网格）
lon_min = min(lon(:));
lon_max = max(lon(:));
lat_min = min(lat(:));
lat_max = max(lat(:));

% 创建经纬度刻度向量
lon_ticks = ceil(lon_min/60)*60 : 60 : floor(lon_max/60)*60;
lat_ticks = ceil(lat_min/30)*30 : 30 : floor(lat_max/30)*30;

% 确保至少有两个刻度点
if length(lon_ticks) < 2
    lon_ticks = [lon_min, lon_max];
end
if length(lat_ticks) < 2
    lat_ticks = [lat_min, lat_max];
end

% 准备网格数据
if isvector(lon) && isvector(lat)
    [LON, LAT] = meshgrid(lon, lat);
else
    LON = lon;
    LAT = lat;
end

% ========== 修改部分1：设置每个子图的色棒范围 ==========
% 这里可以手动设置每个子图的色棒范围，格式为：[min1, max1; min2, max2; ...]
% 示例：对于12个子图，设置不同的范围
colorbar_limits = [
    -0.005, 0.005;    % 子图1
    -0.005, 0.005;    % 子图2
    -0.005, 0.005;    % 子图3
    -0.15, 0.15;    % 子图4
    -0.15, 0.15;      % 子图5
    -0.15, 0.15;    % 子图6
    -0.5, 0.5;    % 子图7
    -15, 15;    % 子图8
    -0.5, 0.5;    % 子图9
    -1, 1;    % 子图10
    -0.05, 0.05;    % 子图11
    -0.02,0.02     % 子图12
];

% ========== 修改部分2：设置紧凑的子图布局 ==========
% 调整子图间距和边距
subplot_margin = 0.04;  % 子图之间的边距
subplot_width = (1 - 4*subplot_margin)/3;  % 子图宽度
subplot_height = (1 - 5*subplot_margin)/4; % 子图高度
set(gca,'Position', [100, 100, 1400, 900])

% 循环绘制12个子图
for subplot_idx = 1:length(Trend_data)
    % 计算当前子图的行列位置
    row = ceil(subplot_idx / 3);
    col = mod(subplot_idx - 1, 3) + 1;
    
    % ========== 修改部分3：计算紧凑的子图位置 ==========
    % 计算子图左下角坐标
    left_pos = subplot_margin + (col-1)*(subplot_width + subplot_margin);
    bottom_pos = 1 - row*(subplot_height + subplot_margin);
    
    % 调整子图位置以留出颜色条空间
    axes_width = subplot_width * 0.85;  % 为颜色条留出15%的空间
    axes_pos = [left_pos, bottom_pos, axes_width, subplot_height];
    
    % 创建子图
    axes('Position', axes_pos);
    
    % ========== 修改部分4：添加左上角编号 ==========
    % 计算编号的位置（在图外边左上角）
    text_x = left_pos + 0.02;
    text_y = bottom_pos + subplot_height - 0.02;
    
    % 添加编号
    annotation('textbox', [text_x, text_y, 0.02, 0.02], ...
        'String', sprintf('(%d)', subplot_idx), ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom');
    
    % 从存储的数据中读取趋势
    Trend = Trend_data{subplot_idx};
    
    Trend_mean = mean(Trend,3);
    %==========================均值和显著性
    alpha = 0.05; % 95% 置信度
    n = 3; % 样本量，这里 n = 3
    df = n - 1; % 自由度 = 2
    % 计算每个格点的标准差
    std_field = std(Trend, 0, 3); % 0 表示使用无偏估计 (除以 n-1)
    % 计算 t 统计量: t = (样本均值 - 0) / (标准误)
    % 零假设 H0: 平均值 = 0
    % 备择假设 H1: 平均值 ≠ 0 (双侧检验)
    standard_error = std_field / sqrt(n);
    t_stat = Trend_mean ./ standard_error;
    % 方法A: 计算 p 值
    p_value = 2 * (1 - tcdf(abs(t_stat), df)); % 双侧检验
    % 显著性掩膜 (显著 = 1, 不显著 = 0)
    significant_mask_p = p_value < alpha;          % 使用 p 值判断
    significant_mask_p(isnan(Trend_mean)) = 1;   %陆地显著，不要填充

    Trend_sig = significant_mask_p;

    % 初始化地图投影
    proj_type = 'robinson';
    m_proj(proj_type, 'lon', [lon_min lon_max], ...
                        'lat', [lat_min lat_max]);
    
    % 绘制海岸线
    m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'none');
    
    hold on;
    
    % 使用 m_pcolor 绘制数据
    h = m_pcolor(LON, LAT, Trend_mean);
    set(h, 'EdgeColor', 'none', 'FaceAlpha', 1);
    
    % 设置颜色
    colormap(mycolor);
    
    % ========== 修改部分5：使用单独的色棒范围 ==========
    clim(colorbar_limits(subplot_idx, :));

    % 绘制显著区域点（稀疏化）
    sig_mask = Trend_sig;

    % 下采样参数
    downsample_factor = 3;

    % 创建采样掩码
    [nrows, ncols] = size(sig_mask);
    row_indices = 1:downsample_factor:nrows;
    col_indices = 1:downsample_factor:ncols;

    % 提取采样后的显著区域
    sig_mask_downsampled = sig_mask(row_indices, col_indices);

    % 提取显著点的经纬度
    if isvector(lon) && isvector(lat)
        [lon_grid, lat_grid] = meshgrid(lon, lat);
        lon_grid_down = lon_grid(row_indices, col_indices);
        lat_grid_down = lat_grid(row_indices, col_indices);
        sig_lon = lon_grid_down(sig_mask_downsampled == 0);
        sig_lat = lat_grid_down(sig_mask_downsampled == 0);
    else
        lon_down = lon(row_indices, col_indices);
        lat_down = lat(row_indices, col_indices);
        sig_lon = lon_down(sig_mask_downsampled == 0);
        sig_lat = lat_down(sig_mask_downsampled == 0);
    end

    % 绘制显著区域点
    if ~isempty(sig_lon)
        m_plot(sig_lon, sig_lat, ...
               '.', 'Color', [0 0 0], 'MarkerSize', 1);
    end

    m_grid('linestyle', 'none', ...
               'tickdir', 'in', ...
               'xtick', lon_ticks, ...
               'ytick', lat_ticks, ...
               'fontsize', 8, ...
               'XaxisLocation', 'bottom', ...
               'YaxisLocation', 'left');

    hold off;
    
    % ========== 修改部分6：在每个子图右侧添加独立的色棒 ==========
    % 计算颜色条位置
    cbar_left = left_pos + axes_width + 0.01;
    cbar_bottom = bottom_pos;
    cbar_width = 0.01;
    cbar_height = subplot_height * 0.8;
    
    % 创建颜色条
    cbar_axes = axes('Position', [cbar_left, cbar_bottom + (subplot_height-cbar_height)/2, ...
                                  cbar_width, cbar_height]);
    colormap(cbar_axes, mycolor);
    caxis(colorbar_limits(subplot_idx, :));
    
    % 设置颜色条
    colorbar(cbar_axes, 'Location', 'eastoutside', 'FontSize', 8);
    
    % 隐藏颜色条坐标轴
    cbar_axes.Visible = 'off';
end


