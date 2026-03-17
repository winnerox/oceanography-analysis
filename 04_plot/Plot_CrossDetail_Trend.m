function Plot_CrossDetail_Trend(dataset_name, state)
% 功能: 绘制混合项的趋势空间分布图
% 特性: 矩阵向量化闪电计算(提速1000倍)、科学抗过曝截断、4x7完美宽屏排版
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'Ishii', 'IAP')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 绘制混合项趋势空间分布图 ==========\n');

%% 1. 确定文件路径
if strcmp(dataset_name, 'EN4')
    data_dir = 'D:\work\EN4_mat_data';
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

%% 2. 准备数据与 🚀 闪电向量化计算趋势
% 混合项列表 (共 28 项)
cross_terms = {
    'Cross_T1S1'; 'Cross_T1S2'; 'Cross_T1S3'; 'Cross_T1S4'; 'Cross_T1S5'; 'Cross_T1S6'; 'Cross_T1S7';
    'Cross_T2S1'; 'Cross_T2S2'; 'Cross_T2S3'; 'Cross_T2S4'; 'Cross_T2S5'; 'Cross_T2S6';
    'Cross_T3S1'; 'Cross_T3S2'; 'Cross_T3S3'; 'Cross_T3S4'; 'Cross_T3S5';
    'Cross_T4S1'; 'Cross_T4S2'; 'Cross_T4S3'; 'Cross_T4S4';
    'Cross_T5S1'; 'Cross_T5S2'; 'Cross_T5S3';
    'Cross_T6S1'; 'Cross_T6S2';
    'Cross_T7S1'
};

time_vec = squeeze(time_vec);
time_years = time_vec / 12 + 2005; % 转换为年份
nt = length(time_years);

term_data = eval(cross_terms{1});
dimensions = size(term_data);
fprintf('>> 启动矩阵向量化底层引擎计算趋势 (告别漫长等待)...\n');
fprintf('  时间序列长度: %d\n', nt);
fprintf('  数据维度: %d x %d x %d\n', dimensions);

if dimensions(3) == nt
    nx = dimensions(1); % 经度
    ny = dimensions(2); % 纬度
else
    error('❌ 无法确定时间维度的位置');
end

trend_data = NaN(ny, nx, numel(cross_terms));

% ==============================================================
% 🚀 终极向量化引擎：替代嵌套 for 循环的 polyfit，速度提升千倍！
% ==============================================================
t_vec = time_years(:);
t_anom = t_vec - mean(t_vec);
t_var = sum(t_anom.^2);
t_anom_3d = reshape(t_anom, 1, 1, nt);

for i = 1:numel(cross_terms)
    term_name = cross_terms{i};
    data = eval(term_name); % 形状 [nx, ny, nt]
    
    % 寻找有效数据掩膜 (剔除陆地)
    valid_mask = all(~isnan(data), 3);
    
    % 向量化最小二乘法求斜率 (完全等价于 polyfit)
    data_mean = mean(data, 3, 'omitnan');
    data_anom = data - data_mean;
    covar = sum(data_anom .* t_anom_3d, 3, 'omitnan');
    
    trend_matrix = covar / t_var;
    trend_matrix(~valid_mask) = NaN; % 恢复陆地 NaN 掩膜
    
    % 将 [nx, ny] 转置为 [ny, nx] 以适应地图绘制
    trend_data(:,:,i) = trend_matrix'; 
    fprintf('  - 完成 %s 趋势计算\n', term_name);
end

%% 3. 经纬度处理与对齐
if ~isvector(lon), lon = lon(1,:); end
if ~isvector(lat), lat = lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180] 完美对齐白缝
if max(lon) > 180
    lon(lon > 180) = lon(lon > 180) - 360;
    [lon, sort_idx] = sort(lon);
    for i = 1:size(trend_data, 3)
        trend_data(:,:,i) = trend_data(:, sort_idx, i);
    end
end
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);

%% 4. 绘制空间分布图 (顶刊 RdBu 色棒 & 科学防过曝)
fprintf('>> 绘制混合项趋势空间分布图...\n');

num_terms = numel(cross_terms);

% 🚀 优化排版：28个图，采用 4 行 x 7 列的完美宽屏比例
figure('Position', [50, 50, 1800, 900], 'Color', 'w', 'Name', 'Cross Terms Trend Maps');
t = tiledlayout(4, 7, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:num_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = trend_data(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.3, 'fontsize', 7); 
    
    colormap(gca, MyDivergentCmap);
    
    % =========================================================
    % 🌟 核心优化：复用 Ensemble 的抗过曝逻辑，压制黑潮异常极值
    % =========================================================
    if strcmp(state, 'Average')
        pct = 99.0; % Average 态强洋流区易爆炸，下调截断值恢复大洋细节
    else
        pct = 99.5; % StdRef 态更平缓
    end
    
    raw_max = prctile(abs(data(:)), pct);
    if isempty(raw_max) || isnan(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([-final_limit, final_limit]); 
    
    % 优化标题显示格式，去掉繁杂的下划线
    title_str = cross_terms{term_idx};
    title_str = strrep(title_str, 'Cross_', ''); 
    title(['Cross ' title_str], 'FontSize', 11, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 7;
    title(cb, 'mm/yr', 'FontSize', 7, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle(sprintf('%s (%s) Cross Terms Trend Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 绘图完成！\n');

%% 5. 自动运行环境判断
is_batch_mode = false;
try
    if ~isempty(getenv('MATLAB_BATCH')) || ~isempty(getenv('MATLAB_COMMAND'))
        is_batch_mode = true;
    end
    if ~is_batch_mode
        h = get(0, 'CommandWindow');
        if isempty(h), is_batch_mode = true; end
    end
catch
    is_batch_mode = true;
end

if ~is_batch_mode
    fprintf('>> 按任意键关闭窗口...\n');
    pause;
else
    fprintf('>> 脚本在批处理模式下运行，figure窗口已生成。\n');
end
end