function Map_Amp_EN4()
% 功能: 绘制 EN4 (Average) 独立项振幅空间分布 (完美 SCI 样式)
% 修复: 
%   1. 【解决过曝】分位数阈值从 95% 提高到 99.5%，还原高振幅区细节。
%   2. 【更换色棒】采用 SCI 专业的平滑“白-黄-红”顺序色盘。
%   3. 【智能刻度】保证色棒刻度整洁。

clear; clc; close all;
fprintf('========== 开始绘制 EN4 振幅分布 (Style Refined) ==========\n');

%% 1. 基础配置
% 你原代码里的路径和文件名
files = {'EN4_TSLA_Terms_1to8_Average.mat'};
input_file = '';
DataDir = 'D:\work\EN4_TSLA_Terms'; % 如果在本地目录请把这个路径改为你存放 .mat 的真实路径，或者把 .mat 放脚本旁边

% 智能加载文件
if exist(files{1}, 'file')
    input_file = files{1};
elseif exist(fullfile(DataDir, files{1}), 'file')
    input_file = fullfile(DataDir, files{1});
else
    error('❌ 找不到数据文件: %s', files{1});
end

fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

%% 2. 坐标对齐与矩阵重排 (防错位)
% 经纬度修正逻辑
if ~exist('lon', 'var'), if exist('Lon','var'), lon=Lon; elseif exist('longitude','var'), lon=longitude; end; end
if ~exist('lat', 'var'), if exist('Lat','var'), lat=Lat; elseif exist('latitude','var'), lat=latitude; end; end
if isvector(lon), [Lon_Grid, Lat_Grid] = meshgrid(lon, lat); else, Lon_Grid=lon; Lat_Grid=lat; lon=Lon_Grid(1,:)'; end
if isvector(lat), lat=Lat_Grid(:,1); end

% 坐标转换 [0, 360] -> [-180, 180]
if max(lon) > 180
    lon_new = lon;
    lon_new(lon_new > 180) = lon_new(lon_new > 180) - 360;
    [lon_sorted, sort_idx] = sort(lon_new);
    [Lon_Grid_New, Lat_Grid_New] = meshgrid(lon_sorted, lat);
else
    Lon_Grid_New = Lon_Grid; Lat_Grid_New = Lat_Grid;
    sort_idx = 1:length(lon);
end

% 处理时间轴 compatibility
if ~exist('time_axis', 'var') && exist('Time_Axis', 'var')
    time_axis = Time_Axis;
end
t_vec = time_axis(:);

%% 3. 核心计算：振幅 (Residual STD)
[~, ~, ntime, max_order] = size(TSLA_AllOrders);
Map_Amp = nan(length(lat), length(lon), 8);

fprintf('>> 极速计算 1-8 阶振幅 (mm)...\n');
for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    % 维度适配
    if size(slice,1) ~= length(lon), slice = permute(slice, [2,1,3]); end
    
    Y = reshape(slice, length(lon)*length(lat), ntime)';
    valid_idx = ~isnan(Y(1,:));
    Y_valid = Y(:, valid_idx);
    
    % 去趋势最小二乘
    X_reg = [t_vec - mean(t_vec), ones(ntime, 1)];
    B = X_reg \ Y_valid;
    Res = Y_valid - (X_reg * B);
    
    amp_valid = std(Res, 0, 1);
    amp_1d = nan(1, length(lon)*length(lat));
    amp_1d(valid_idx) = amp_valid;
    temp_map = reshape(amp_1d, length(lon), length(lat))';
    
    % 应用经度重排
    Map_Amp(:,:,n) = temp_map(:, sort_idx);
end

% 单位统一为 mm
if max(Map_Amp(:)) < 10, Map_Amp = Map_Amp * 1000; end

%% 4. ==========================================
%  绘图配置：复刻 Style 2 样式
%  ==========================================
figure('Position', [100, 50, 1100, 1300], 'Color', 'w'); 
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% ======= 【色棒修改】专业 SCI 暖色系顺序色盘 (白-黄-红) =======
MyRedPoints = [
    1,     1,     1;      % 0.0: 纯白 (极弱区域)
    1,     1,     0.7;    % 0.1: 暖白
    1,     0.9,   0.5;    % 0.4: 淡黄
    0.95,  0.6,   0.2;    % 0.8: 橙红
    0.8,   0.1,   0.1;    % 1.0: 红 (强振幅区)
    0.5,   0,     0       % >1 : 深红
];
xx = [0, 0.1, 0.4, 0.7, 0.9, 1]; % 节点位置，保持白色和黄色区较宽，红色区紧凑
yy = linspace(0,1,256);
MySequentialCmap = [interp1(xx,MyRedPoints(:,1),yy,'pchip')', ...
                    interp1(xx,MyRedPoints(:,2),yy,'pchip')', ...
                    interp1(xx,MyRedPoints(:,3),yy,'pchip')'];

for term_idx = 1:8
    ax = nexttile;
    
    % 地图投影
    m_proj('robinson', 'lon', [-180 180], 'lat', [-80 80]); % 修正Robinson南极过窄问题
    data = Map_Amp(:,:,term_idx);
    
    % 绘制填色图
    m_pcolor(Lon_Grid_New, Lat_Grid_New, data); shading flat; hold on;
    
    % 绘制海岸线和网格
    m_coast('patch', [.96 .96 .96], 'edgecolor', [.5 .5 .5], 'linewidth', 0.5);
    
    % ======= 【网格刻度修改】强制显色、细化经纬度 =======
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.75 .75 .75], ...
           'linewidth', 0.5, 'fontsize', 8, 'XaxisLocation', 'bottom', 'YaxisLocation', 'left'); 
    
    % 应用新的色盘
    colormap(ax, MySequentialCmap);
    
    % ======= 【解决过曝】高阶分位数取值策略 =======
    % 从 95% 提高到 99.5%，极大化还原细节，防止死红
    raw_max = prctile(data(:), 99.5); 
    if isnan(raw_max) || raw_max < 1e-12, raw_max = 1; end
    
    % 智能刻度取整
    exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
    if fraction <= 1.5, nice_b = 1.5; elseif fraction <= 3, nice_b = 3; elseif fraction <= 5, nice_b = 5; else, nice_b = 10; end
    final_l = nice_b * 10^exponent;
    
    % 设置范围
    caxis([0, final_l]); 
    
    title(sprintf('Term %d Amplitude', term_idx), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    
    % 色棒样式调整
    cb = colorbar; cb.FontSize = 8;
    title(cb, 'mm', 'FontSize', 8, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    cb.Ticks = linspace(0, final_l, 5); % 保持5个整洁刻度
    try cb.Ruler.Exponent = 0; catch; end
end

% 大标题
sgtitle('\fontname{Microsoft YaHei}EN4 (Standard Ref.) 各阶泰勒展开项振幅空间分布', 'FontSize', 18, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto'); 

fprintf('\n>> 绘制完成！过曝已解决，色棒已改为 SCI 专业 Style 2。\n');
end