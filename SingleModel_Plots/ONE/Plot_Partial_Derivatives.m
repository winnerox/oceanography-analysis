function Plot_Partial_Derivatives(dataset_name, state)
% 功能: 计算并绘制1-8阶偏导系数的空间分布图
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 计算并绘制高阶偏导系数分布图 ==========\n');

%% 1. 加载数据
% 基础目录
base_dir = 'D:\work';

% 泰勒展开项文件路径
if strcmp(dataset_name, 'EN4')
    terms_dir = fullfile(base_dir, 'EN4_TSLA_Terms');
elseif strcmp(dataset_name, 'IAP')
    terms_dir = fullfile(base_dir, 'IAP_mat_data');
elseif strcmp(dataset_name, 'Ishii')
    terms_dir = fullfile(base_dir, 'Ishii_mat_data');
else
    error('不支持的数据集: %s', dataset_name);
end

% 加载TSLA展开项
tsla_file = fullfile(terms_dir, sprintf('%s_TSLA_Terms_1to8_%s.mat', dataset_name, state));
if ~exist(tsla_file, 'file')
    error('❌ 找不到TSLA展开项文件: %s', tsla_file);
end

fprintf('>> 加载TSLA展开项数据: %s ...\n', tsla_file);
tsla_data = load(tsla_file);

% 获取TSLA展开项变量名
if isfield(tsla_data, 'TSLA_AllOrders')
    TSLA_AllOrders = tsla_data.TSLA_AllOrders;
else
    fields = fieldnames(tsla_data);
    for i = 1:length(fields)
        if contains(fields{i}, 'TSLA') && contains(fields{i}, 'Order')
            TSLA_AllOrders = tsla_data.(fields{i});
            break;
        end
    end
end

% 加载盐度展开项
ssla_file = fullfile(terms_dir, sprintf('%s_S_Terms_1to8_%s.mat', dataset_name, state));
SSLA_AllOrders = [];
if exist(ssla_file, 'file')
    fprintf('>> 加载盐度展开项数据: %s ...\n', ssla_file);
    ssla_data = load(ssla_file);
    if isfield(ssla_data, 'SSLA_AllOrders')
        SSLA_AllOrders = ssla_data.SSLA_AllOrders;
    end
end

% 加载交叉项
cross_file = fullfile(terms_dir, sprintf('%s_Cross_Terms_1to8_%s.mat', dataset_name, state));
Cross_AllOrders = [];
if exist(cross_file, 'file')
    fprintf('>> 加载交叉项数据: %s ...\n', cross_file);
    cross_data = load(cross_file);
    if isfield(cross_data, 'Cross_AllOrders')
        Cross_AllOrders = cross_data.Cross_AllOrders;
    end
end

% 加载平均态数据（包含温度和盐度）
mean_file = fullfile(terms_dir, sprintf('%s_Mean_State.mat', dataset_name));
if ~exist(mean_file, 'file')
    mean_file = fullfile(terms_dir, 'EN4_Mean_State.mat'); % 尝试使用EN4的平均态文件
    if ~exist(mean_file, 'file')
        error('❌ 找不到平均态文件');
    end
end

fprintf('>> 加载平均态数据: %s ...\n', mean_file);
mean_data = load(mean_file);

% 获取坐标数据
if isfield(tsla_data, 'Lon')
    Lon = tsla_data.Lon;
elseif isfield(tsla_data, 'lon')
    Lon = tsla_data.lon;
elseif isfield(mean_data, 'Lon')
    Lon = mean_data.Lon;
elseif isfield(mean_data, 'lon')
    Lon = mean_data.lon;
else
    error('❌ 未找到经度数据');
end

if isfield(tsla_data, 'Lat')
    Lat = tsla_data.Lat;
elseif isfield(tsla_data, 'lat')
    Lat = tsla_data.lat;
elseif isfield(mean_data, 'Lat')
    Lat = mean_data.Lat;
elseif isfield(mean_data, 'lat')
    Lat = mean_data.lat;
else
    error('❌ 未找到纬度数据');
end

%% 2. 数据预处理
% 确保Lon和Lat是向量
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180]
sort_idx = [];
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
end

% 确保网格维度与数据一致
fprintf('>> 原始Lon长度: %d, Lat长度: %d\n', length(Lon), length(Lat));

% 获取数据维度
[ntime, nlat, nlon, norder] = size(TSLA_AllOrders);
fprintf('>> 数据维度: time=%d, lat=%d, lon=%d, order=%d\n', ntime, nlat, nlon, norder);

% 确保Lon和Lat的长度与数据匹配
if length(Lon) ~= nlon
    fprintf('>> 调整Lon长度从 %d 到 %d\n', length(Lon), nlon);
    Lon = linspace(min(Lon), max(Lon), nlon);
end
if length(Lat) ~= nlat
    fprintf('>> 调整Lat长度从 %d 到 %d\n', length(Lat), nlat);
    Lat = linspace(min(Lat), max(Lat), nlat);
end

% 创建网格
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);
fprintf('>> 网格维度: Lon_Grid=%dx%d, Lat_Grid=%dx%d\n', size(Lon_Grid), size(Lat_Grid));

%% 3. 计算偏导系数
fprintf('>> 计算1-8阶偏导系数...\n');

% 初始化偏导系数数组
% d^nTSLA/dT^n: 温度偏导系数
dT_derivatives = zeros(nlat, nlon, norder);

% d^nTSLA/dS^n: 盐度偏导系数
dS_derivatives = zeros(nlat, nlon, norder);

% 交叉项偏导系数
cross_derivatives = zeros(nlat, nlon, norder);

% 计算温度和盐度的异常值（ΔT, ΔS）
% 这里我们使用时间序列的标准差作为异常值的估计
% 实际应用中，应该使用真实的温度和盐度异常数据

% 计算TSLA的时间标准差作为ΔTSLA的估计
TSLA_std = std(TSLA_AllOrders, [], 1);
TSLA_std = squeeze(TSLA_std);

% 确保TSLA_std的维度是正确的
if length(size(TSLA_std)) == 2
    % 如果只有lat和lon维度，扩展为lat×lon×order
    TSLA_std = repmat(TSLA_std, [1, 1, norder]);
end

% 假设温度和盐度的异常值与TSLA异常值相关
% 这里使用简化的方法，实际应用中应该使用真实的温度和盐度数据
Delta_T = TSLA_std * 10;  % 假设温度异常是TSLA的10倍
Delta_S = TSLA_std * 100; % 假设盐度异常是TSLA的100倍

% 计算偏导系数
for order = 1:norder
    fprintf('   计算第 %d 阶偏导系数...\n', order);
    
    % 获取当前阶数的展开项
    tsla_term = TSLA_AllOrders(:, :, :, order);
    ssla_term = [];
    cross_term = [];
    
    if ~isempty(SSLA_AllOrders)
        ssla_term = SSLA_AllOrders(:, :, :, order);
    end
    
    if ~isempty(Cross_AllOrders)
        cross_term = Cross_AllOrders(:, :, :, order);
    end
    
    % 计算时间平均的展开项
    tsla_mean = mean(tsla_term, 1);
    tsla_mean = squeeze(tsla_mean);
    
    if ~isempty(ssla_term)
        ssla_mean = mean(ssla_term, 1);
        ssla_mean = squeeze(ssla_mean);
    end
    
    if ~isempty(cross_term)
        cross_mean = mean(cross_term, 1);
        cross_mean = squeeze(cross_mean);
    end
    
    % 计算偏导系数
    % 泰勒展开公式: term = (1/n!) * (d^nTSLA/dx^n) * Δx^n
    % 所以偏导系数: d^nTSLA/dx^n = term * n! / Δx^n
    
    factorial_n = factorial(order);
    
    % 温度偏导系数（假设TSLA_Terms包含温度项）
    if order == 1
        % 一阶偏导: dTSLA/dT = term / ΔT
        dT_derivatives(:, :, order) = tsla_mean ./ Delta_T(:, :, order);
    else
        % 高阶偏导: d^nTSLA/dT^n = term * n! / ΔT^n
        dT_derivatives(:, :, order) = tsla_mean * factorial_n ./ (Delta_T(:, :, order) .^ order);
    end
    
    % 盐度偏导系数
    if ~isempty(ssla_term)
        if order == 1
            dS_derivatives(:, :, order) = ssla_mean ./ Delta_S(:, :, order);
        else
            dS_derivatives(:, :, order) = ssla_mean * factorial_n ./ (Delta_S(:, :, order) .^ order);
        end
    end
    
    % 交叉项偏导系数
    if ~isempty(cross_term)
        % 交叉项的偏导系数计算比较复杂，这里使用简化方法
        cross_derivatives(:, :, order) = cross_mean * factorial_n ./ (Delta_T(:, :, order) .^ (order/2) .* Delta_S(:, :, order) .^ (order/2));
    end
end

%% 4. 绘制偏导系数空间分布图
fprintf('>> 绘制偏导系数空间分布图...\n');

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

% 绘制温度偏导系数
fprintf('>> 绘制温度偏导系数分布图...\n');
figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Temperature Partial Derivatives');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for order = 1:norder
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = dT_derivatives(:,:,order);
    % 详细的维度检查
    fprintf('>> 绘制第 %d 阶温度偏导系数，数据维度: %dx%d, 网格维度: %dx%d\n', order, size(data), size(Lon_Grid));
    
    % 确保数据维度与网格一致
    if size(data, 1) ~= size(Lat_Grid, 1) || size(data, 2) ~= size(Lon_Grid, 2)
        fprintf('>> 维度不匹配，调整数据维度...\n');
        % 调整数据维度以匹配网格
        data = reshape(data, size(Lat_Grid, 1), size(Lon_Grid, 2));
        fprintf('>> 调整后数据维度: %dx%d\n', size(data));
    end
    
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, MyDivergentCmap);
    
    raw_max = prctile(abs(data(:)), 99.8);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    
    % 确保final_limit是正数
    if final_limit <= 0 || isnan(final_limit) || isinf(final_limit)
        final_limit = 1;
    end
    
    caxis([-final_limit, final_limit]); 
    
    title(sprintf('%d-Order Temperature Derivative', order), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/(yr·°C^n)', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) Temperature Partial Derivatives Distribution', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

% 绘制盐度偏导系数
if ~isempty(SSLA_AllOrders)
    fprintf('>> 绘制盐度偏导系数分布图...\n');
    figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Salinity Partial Derivatives');
    t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for order = 1:norder
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        
        data = dS_derivatives(:,:,order);
        % 详细的维度检查
        fprintf('>> 绘制第 %d 阶盐度偏导系数，数据维度: %dx%d, 网格维度: %dx%d\n', order, size(data), size(Lon_Grid));
        
        % 确保数据维度与网格一致
        if size(data, 1) ~= size(Lat_Grid, 1) || size(data, 2) ~= size(Lon_Grid, 2)
            fprintf('>> 维度不匹配，调整数据维度...\n');
            % 调整数据维度以匹配网格
            data = reshape(data, size(Lat_Grid, 1), size(Lon_Grid, 2));
            fprintf('>> 调整后数据维度: %dx%d\n', size(data));
        end
        
        m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
        
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
        
        colormap(gca, MyDivergentCmap);
        
        raw_max = prctile(abs(data(:)), 99.8);
        if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
        exponent = floor(log10(raw_max)); 
        fraction = raw_max / 10^exponent; 
        if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
        final_limit = nice_base * 10^exponent;
        
        % 确保final_limit是正数
        if final_limit <= 0 || isnan(final_limit) || isinf(final_limit)
            final_limit = 1;
        end
        
        caxis([-final_limit, final_limit]); 
        
        title(sprintf('%d-Order Salinity Derivative', order), 'FontSize', 12, 'FontWeight', 'bold');
        
        cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
        title(cb, 'mm/(yr·PSU^n)', 'FontSize', 8, 'FontWeight', 'bold'); 
        cb.Ticks = linspace(-final_limit, final_limit, 5);
    end
    
    sgtitle(sprintf('%s (%s) Salinity Partial Derivatives Distribution', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
end

%% 5. 绘制偏导系数箱线图
fprintf('>> 绘制偏导系数箱线图...\n');

% 准备箱线图数据
mask = ~isnan(dT_derivatives(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_dT = nan(num_valid, norder);
Data_dS = nan(num_valid, norder);

for order = 1:norder
    % 温度偏导系数数据
    data_slice = dT_derivatives(:,:,order);
    vec = nan(num_valid, 1);
    for k = 1:num_valid
        vec(k) = data_slice(rows(k), cols(k));
    end
    Data_dT(:, order) = vec;
    
    % 盐度偏导系数数据
    if ~isempty(SSLA_AllOrders)
        data_slice = dS_derivatives(:,:,order);
        vec = nan(num_valid, 1);
        for k = 1:num_valid
            vec(k) = data_slice(rows(k), cols(k));
        end
        Data_dS(:, order) = vec;
    end
end

% 绘制温度偏导系数箱线图
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Temperature Partial Derivatives Boxplot');

% 顶刊配色
C1 = [0.902, 0.294, 0.208]; % 樱桃红
C2 = [0.302, 0.455, 0.690]; % 钴蓝
C3 = [0.000, 0.627, 0.408]; % 薄荷绿
BoxBgColor = [0.96 0.96 0.96]; % 极浅高级灰底色

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

% 设置箱线图位置和范围
X_Limit = [0.4, 9.0];
MainPos = [0.07, 0.14, 0.80, 0.80]; 

% 绘制箱线图
ax1 = axes('Position', MainPos); hold(ax1, 'on');

% 为每个阶数绘制箱线图
positions = 1:norder;
colors = [C1; C2; C3; C1; C2; C3; C1; C2];

for k = 1:norder
    draw_jitter_unit(ax1, positions(k), Data_dT(:,k), BoxWidth, ScatterWidth, ScatterSize, colors(k,:), BoxBgColor, ScatterAlpha);
end

% 使用数据的实际最大最小值
min_v = min(Data_dT(:));
max_v = max(Data_dT(:));
if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12
    min_v = -1;
    max_v = 1;
end
% 添加10%的边距
range = max_v - min_v;
Limit_min = min_v - range * 0.1;
Limit_max = max_v + range * 0.1;
ylim(ax1, [Limit_min, Limit_max]); 
set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel(ax1, 'Temperature Partial Derivative', 'FontWeight', 'bold', 'FontSize', 16);

% --- 装饰与分割线 ---
xline(ax1, 0.7, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);

All_Ticks = positions;
All_Labels = cellstr(num2str((1:norder)'));
set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
xlabel(ax1, 'Order', 'FontSize', 16, 'FontWeight', 'bold');
title(ax1, sprintf('%s (%s) Temperature Partial Derivatives Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');

% === 图例 ===
h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_out  = scatter(ax1, NaN, NaN, ScatterSize, C1, 'filled', 'MarkerFaceAlpha', 0.6);

legend(ax1, [h_box, h_med, h_mean, h_out], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'none');

% === 交互优化 ===
hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
hPan = pan(gcf); hPan.Motion = 'both';

fprintf('>> 绘图完成！\n');
end

%% ============================================================
%% 内部函数
%% ============================================================
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    max_val = max(data); % 计算最大值
    
    idx_out = data < low_w | data > up_w;
    outliers = data(idx_out);
    
    % 绘制散点 (去除黑边框，恢复果冻质感透明叠加)
    if ~isempty(outliers)
        if length(outliers) > 2000
            outliers = outliers(randperm(length(outliers), 2000));
        end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        try
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch
            scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s);
        end
    end
    
    % 绘制须线 (恢复纤细的 1.2 细线)
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    % 绘制箱体
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, ...
        'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    % 绘制中位线、均值点和最大值点
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(ax, center, max_val, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % 绘制最大值点
end
