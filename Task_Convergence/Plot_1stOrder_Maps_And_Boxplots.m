function Plot_1stOrder_Maps_And_Boxplots()
% 功能：基于用户原始代码风格 (draw_jitter_unit & m_map) 
% 绘制三套数据 全部八阶(HSLA) 的箱线图与空间分布图
% 修复：处理了 0~360 经度到 -180~180 的转换与矩阵重排，防止地图错位偏移

clear; clc; close all;
fprintf('========== 开始绘制全部八阶对比图 (Style Matched) ==========\n');

% --- 1. 基础设置与数据加载 ---
base_dir = 'D:\work';
trend_dir = fullfile(base_dir, 'Task_Convergence', 'Trend_Results');

data_names = {'EN4', 'IAP', 'Ishii'};
state = 'Average'; % 如果画标准态改这里: 'StdRef'
max_order = 8;     % 绘制全部8阶

trend_1d = cell(1,3); amp_1d = cell(1,3);
t_maps = cell(max_order, 3);   a_maps = cell(max_order, 3);   s_maps = cell(max_order, 3);
lon_grids = cell(1,3); lat_grids = cell(1,3);

for i = 1:3
    ds_name = data_names{i};
    file_path = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', ds_name, state));
    if ~exist(file_path, 'file'), error('找不到文件: %s', file_path); end
    
    fprintf('>> 读取数据: %s\n', ds_name);
    d = load(file_path);
    
    Lon = d.Lon; Lat = d.Lat;
    
    % --- 规范化一维经纬度向量 ---
    if isvector(Lon)
        lon_vec = Lon(:); lat_vec = Lat(:);
    else
        lon_vec = Lon(1,:)'; lat_vec = Lat(:,1);
    end
    
    % =======================================================
    % 🚀 核心修复：经纬度 [0, 360] 转换到 [-180, 180] 并重排矩阵
    % =======================================================
    if max(lon_vec) > 180
        lon_vec(lon_vec > 180) = lon_vec(lon_vec > 180) - 360;
    end
    [lon_vec, sort_idx] = sort(lon_vec); % 获取正确的经度排序索引
    [Lon_Grid, Lat_Grid] = meshgrid(lon_vec, lat_vec);
    
    % 从Amplitude_Results文件夹加载振幅数据
    fprintf('   - 正在加载 %s 振幅数据...\n', ds_name);
    amp_file = fullfile('D:\work\Task_Convergence\Amplitude_Results', sprintf('%s_%s_Amplitude.mat', ds_name, state));
    amp_data = load(amp_file);
    
    % 存储 1D 数组(去NaN)，供箱线图使用
    trend_1d{i} = cell(1, max_order);
    amp_1d{i} = cell(1, max_order);
    
    % 对每个阶数加载数据
    for order = 1:max_order
        % 提取当前阶数的趋势数据 (优先使用总海平面异常 HSLA)
        if isfield(d, 'trend_HSLA') && ~isempty(d.trend_HSLA)
            t_map = d.trend_HSLA(:,:,order);
            s_map = d.sig_HSLA(:,:,order);
        else
            t_map = d.trend_TSLA(:,:,order);
            s_map = d.sig_TSLA(:,:,order);
        end
        
        % 提取当前阶数的振幅数据
        if isfield(amp_data, 'amp_HSLA') && ~isempty(amp_data.amp_HSLA)
            a_map = amp_data.amp_HSLA(:,:,order);
        else
            a_map = amp_data.amp_TSLA(:,:,order);
        end
        
        % 统一矩阵尺寸：确保地图尺寸为 (nlat x nlon)
        if size(t_map, 1) ~= length(lat_vec)
            t_map = t_map';
            s_map = s_map';
        end
        if size(a_map, 1) ~= length(lat_vec)
            a_map = a_map';
        end
        
        % 应用排序索引，将数据跟着经度一起移动，完美居中！
        t_map = t_map(:, sort_idx);
        s_map = s_map(:, sort_idx);
        a_map = a_map(:, sort_idx);
        % =======================================================
        
        % 存储 2D 地图 
        t_maps{order, i} = t_map; a_maps{order, i} = a_map; s_maps{order, i} = s_map;
        
        % 存储 1D 数组(去NaN)，供箱线图使用
        trend_1d{i}{order} = t_map(~isnan(t_map));
        amp_1d{i}{order} = a_map(~isnan(a_map));
    end
    
    lon_grids{i} = Lon_Grid; lat_grids{i} = Lat_Grid;
end

%% ==========================================================
%% Figure 1: 趋势箱线图 (8阶，2x4布局)
%% ==========================================================
fprintf('\n>> 绘制趋势箱线图...\n');
fig1 = figure('Position', [100, 200, 1600, 800], 'Color', 'w', 'Name', 'Trend Boxplots');
colors = {[0.85 0.325 0.098], [0 0.447 0.741], [0.466 0.674 0.188]}; % 红、蓝、绿
BoxBgColor = [0.95 0.95 1];

% 创建2x4布局，显示8阶趋势
for order = 1:8
    ax = subplot(2, 4, order); hold(ax, 'on');
    
    for i = 1:3
        draw_jitter_unit(ax, i, trend_1d{i}{order}, 0.55, 0.4, 12, colors{i}, BoxBgColor, 0.4);
    end
    
    yline(ax, 0, '--k', 'LineWidth', 1.5);
    set(ax, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(ax, 'XTick', 1:3, 'XTickLabel', data_names);
    ylabel(ax, 'Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 12);
    title(ax, sprintf('%dth Order', order), 'FontSize', 14, 'FontWeight', 'bold');
    
    % 计算y轴范围
    max_t = 0;
    for i = 1:3
        if ~isempty(trend_1d{i}{order})
            current_max = prctile(abs(trend_1d{i}{order}), 99.5);
            if current_max > max_t
                max_t = current_max;
            end
        end
    end
    
    if max_t > 0
        ylim(ax, [-max_t*1.2, max_t*1.2]);
    else
        ylim(ax, [-1, 1]);
    end
end



%% ==========================================================
%% Figure 2: 振幅箱线图 (8阶，2x4布局)
%% ==========================================================
fprintf('\n>> 绘制振幅箱线图...\n');
fig2 = figure('Position', [100, 200, 1600, 800], 'Color', 'w', 'Name', 'Amplitude Boxplots');

% 创建2x4布局，显示8阶振幅
for order = 1:8
    ax = subplot(2, 4, order); hold(ax, 'on');
    
    for i = 1:3
        draw_jitter_unit(ax, i, amp_1d{i}{order}, 0.55, 0.4, 12, colors{i}, BoxBgColor, 0.4);
    end
    
    yline(ax, 0, '-k', 'LineWidth', 1.5);
    set(ax, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(ax, 'XTick', 1:3, 'XTickLabel', data_names);
    ylabel(ax, 'Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 12);
    title(ax, sprintf('%dth Order', order), 'FontSize', 14, 'FontWeight', 'bold');
    
    % 计算y轴范围
    max_a = 0;
    for i = 1:3
        if ~isempty(amp_1d{i}{order})
            current_max = prctile(amp_1d{i}{order}, 99.5);
            if current_max > max_a
                max_a = current_max;
            end
        end
    end
    
    if max_a > 0
        ylim(ax, [0, max_a*1.2]);
    else
        ylim(ax, [0, 1]);
    end
end




fprintf('\n>> 完美！全部出图完毕！\n');
end

%% ==========================================================
%% 内部函数库 (完全继承你原始的画图逻辑)
%% ==========================================================

function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    max_val = max(data); % 计算最大值
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    if ~isempty(outliers)
        if length(outliers) > 2000, outliers = outliers(randperm(length(outliers), 2000)); end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        try
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s);
        end
    end
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(ax, center, max_val, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % 绘制最大值点
end

function amp_map = calc_amp_fast(ds_path, ds_name, state, time_vec)
    % 从Amplitude_Results文件夹加载振幅数据
    amp_dir = 'D:\work\Task_Convergence\Amplitude_Results';
    amp_file = fullfile(amp_dir, sprintf('%s_%s_Amplitude.mat', ds_name, state));
    
    if ~exist(amp_file, 'file'), amp_map = nan(360,173); return; end
    
    d = load(amp_file);
    
    % 提取第一阶振幅数据
    if isfield(d, 'amp_HSLA') && ~isempty(d.amp_HSLA)
        amp_map = d.amp_HSLA(:,:,1);
    else
        amp_map = d.amp_TSLA(:,:,1);
    end
end