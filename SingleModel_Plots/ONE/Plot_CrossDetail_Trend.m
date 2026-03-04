function Plot_CrossDetail_Trend(dataset_name, state)
% 功能: 绘制混合项的趋势空间分布图
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 绘制混合项趋势空间分布图 ==========\n');

%% 1. 确定文件路径
if strcmp(dataset_name, 'EN4')
    data_dir = 'D:\work\EN4_TSLA_Terms';
elseif strcmp(dataset_name, 'Ishii')
    data_dir = 'D:\work\Ishii_mat_data';
elseif strcmp(dataset_name, 'IAP')
    data_dir = 'D:\work\IAP_mat_data';
else
    error('❌ 不支持的数据集名称: %s', dataset_name);
end

% 确定文件后缀
if strcmp(state, 'Average')
    state_suffix = 'Avg';
elseif strcmp(state, 'StdRef')
    state_suffix = 'Std';
else
    error('❌ 不支持的状态: %s', state);
end

crossdetail_file = fullfile(data_dir, sprintf('%s_CrossDetail_%s.mat', dataset_name, state_suffix));

if ~exist(crossdetail_file, 'file')
    error('❌ 找不到CrossDetail文件: %s', crossdetail_file);
end

fprintf('>> 加载CrossDetail数据: %s ...\n', crossdetail_file);
load(crossdetail_file);

%% 2. 准备数据
% 混合项列表
cross_terms = {
    'Cross_T1S1'; 'Cross_T1S2'; 'Cross_T1S3'; 'Cross_T1S4'; 'Cross_T1S5'; 'Cross_T1S6'; 'Cross_T1S7';
    'Cross_T2S1'; 'Cross_T2S2'; 'Cross_T2S3'; 'Cross_T2S4'; 'Cross_T2S5'; 'Cross_T2S6';
    'Cross_T3S1'; 'Cross_T3S2'; 'Cross_T3S3'; 'Cross_T3S4'; 'Cross_T3S5';
    'Cross_T4S1'; 'Cross_T4S2'; 'Cross_T4S3'; 'Cross_T4S4';
    'Cross_T5S1'; 'Cross_T5S2'; 'Cross_T5S3';
    'Cross_T6S1'; 'Cross_T6S2';
    'Cross_T7S1'
};

% 计算趋势
time_vec = squeeze(time_vec);
time_years = time_vec / 12 + 2005; % 转换为年份
nt = length(time_years);

% 初始化趋势数据结构
term_data = eval(cross_terms{1});
dimensions = size(term_data);
fprintf('>> 计算混合项趋势...\n');
fprintf('  时间序列长度: %d\n', nt);
fprintf('  数据维度: %d x %d x %d\n', dimensions);

% 确定正确的维度顺序
if dimensions(3) == nt
    % 时间维度在第三维
    time_dim = 3;
    lon_dim = 1; % 经度在第一维
    lat_dim = 2; % 纬度在第二维
    nx = dimensions(lon_dim);
    ny = dimensions(lat_dim);
else
    error('❌ 无法确定时间维度的位置');
end

trend_data = NaN(ny, nx, numel(cross_terms));

for i = 1:numel(cross_terms)
    term_name = cross_terms{i};
    data = eval(term_name);
    
    % 计算每个网格点的趋势
    for y = 1:ny  % 纬度
        for x = 1:nx  % 经度
            time_series = squeeze(data(x, y, :)); % 数据结构: (lon, lat, time)
            ts_length = length(time_series);
            if ts_length == nt && ~all(isnan(time_series))
                p = polyfit(time_years, time_series, 1);
                trend_data(y, x, i) = p(1); % 斜率即为趋势
            elseif ts_length ~= nt
                fprintf('  - 警告: %s 在 (lat=%d, lon=%d) 处时间序列长度不匹配: %d vs %d\n', term_name, y, x, ts_length, nt);
            end
        end
    end
    fprintf('  - 完成 %s 趋势计算\n', term_name);
end

%% 3. 经纬度处理
% 确保Lon和Lat是向量
if ~isvector(lon), lon = lon(1,:); end
if ~isvector(lat), lat = lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180] 完美对齐
if max(lon) > 180
    lon(lon > 180) = lon(lon > 180) - 360;
    [lon, sort_idx] = sort(lon);
    for i = 1:size(trend_data, 3)
        trend_data(:,:,i) = trend_data(:, sort_idx, i);
    end
end
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);

%% 4. 绘制空间分布图 (顶刊 RdBu 色棒)
fprintf('>> 绘制混合项趋势空间分布图...\n');

% 确定显示的布局
num_terms = numel(cross_terms);
rows = ceil(sqrt(num_terms));
cols = ceil(num_terms / rows);

figure('Position', [100, 50, 1200, 800], 'Color', 'w', 'Name', 'Cross Terms Trend Maps');
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

for term_idx = 1:num_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = trend_data(:,:,term_idx);
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
    caxis([-final_limit, final_limit]); 
    
    % 确保标题字符串格式正确
    title_str = cross_terms{term_idx};
    title_str = strrep(title_str, '_', ''); % 移除下划线
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) Cross Terms Trend Distribution', dataset_name, state), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 绘图完成！\n');

% 检查是否在交互式环境中运行
% 在批处理模式下，MATLAB会设置特定的环境变量
is_batch_mode = false;
try
    % 检查是否存在命令行参数
    if ~isempty(getenv('MATLAB_BATCH')) || ~isempty(getenv('MATLAB_COMMAND'))
        is_batch_mode = true;
    end
    % 尝试检查命令窗口是否存在
    if ~is_batch_mode
        h = get(0, 'CommandWindow');
        if isempty(h)
            is_batch_mode = true;
        end
    end
catch
    % 如果出现错误，默认认为是批处理模式
    is_batch_mode = true;
end

if ~is_batch_mode
    fprintf('>> 按任意键关闭窗口...\n');
    pause;
else
    fprintf('>> 脚本在批处理模式下运行，figure窗口已保存。\n');
end
end
