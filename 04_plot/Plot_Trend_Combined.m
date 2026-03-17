function Plot_Trend_Combined(dataset_name, state, target_comp)
% 功能: 整合绘制趋势空间分布图和箱线图 (顶刊标准：不显著区域打点遮盖)
% 特性: 动态自适应阶数，完美兼容 9阶(Exact+8阶Taylor)、8阶(TSLA/HSLA)、单阶 数据
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')
%   target_comp: 要绘制的变量名称 ('TSLA', 'HSLA', 'SSLA') [可选，默认 'SSLA']

clc; close all;
fprintf('========== 整合绘制趋势空间分布图和箱线图 ==========\n');

if nargin < 3
    target_comp = 'SSLA'; % 默认优先尝试画 SSLA 的 9 层结构
end

%% 1. 加载趋势结果文件
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', dataset_name, state));
if ~exist(trend_file, 'file')
    error('❌ 找不到趋势结果文件: %s', trend_file);
end
fprintf('>> 加载趋势数据: %s ...\n', trend_file);
load(trend_file);

%% 2. 准备数据 (🚨 严格的变量隔离与自动合成机制)
fprintf('>> 🎯 当前请求绘制的目标变量: 【%s】\n', target_comp);

if strcmpi(target_comp, 'SSLA')
    if ~exist('trend_SSLA', 'var') || isempty(trend_SSLA)
        error('❌ 请求了 SSLA，但文件中未找到 trend_SSLA 变量！');
    end
    trend_data = trend_SSLA;
    comp_name = 'SSLA (Total Steric)';
    if exist('sig_SSLA', 'var'), sig_data = sig_SSLA; else, sig_data = zeros(size(trend_data)); end
    
    % 🚀 如果同时存在 TSLA 和 HSLA，自动启动 9 层合成引擎
    if exist('trend_TSLA', 'var') && exist('trend_HSLA', 'var')
        fprintf('>> 💡 正在自动合成 1-8 阶泰勒展开的总比容 (TSLA_n + HSLA_n)...\n');
        taylor_ssla = trend_TSLA + trend_HSLA;
        if exist('trend_Cross', 'var')
            taylor_ssla = taylor_ssla + trend_Cross;
        end
        taylor_sig = sig_TSLA; % 借用支配项 TSLA 的显著性作为整体掩膜
        
        % 将 1 层的 Exact 和 8 层的 Taylor 拼成 9 层无缝巨型矩阵！
        trend_data = cat(3, trend_data, taylor_ssla);
        sig_data = cat(3, sig_data, taylor_sig);
        fprintf('>> 🎊 成功生成 9 层总比容对比矩阵 (1层 Exact + 8层 Taylor)！\n');
    else
        fprintf('>> ⚠️ 缺少 TSLA 或 HSLA，仅绘制单层 Exact SSLA。\n');
    end

elseif strcmpi(target_comp, 'TSLA')
    if ~exist('trend_TSLA', 'var') || isempty(trend_TSLA)
        error('❌ 请求了 TSLA，但文件中未找到 trend_TSLA 变量！');
    end
    trend_data = trend_TSLA;
    comp_name = 'TSLA (Thermosteric)';
    if exist('sig_TSLA', 'var'), sig_data = sig_TSLA; else, sig_data = zeros(size(trend_data)); end
    fprintf('>> ✅ 成功提取并锁定纯温度 TSLA 的趋势数据。\n');

elseif strcmpi(target_comp, 'HSLA')
    if ~exist('trend_HSLA', 'var') || isempty(trend_HSLA)
        error('❌ 请求了 HSLA，但文件中未找到 trend_HSLA 变量！');
    end
    trend_data = trend_HSLA;
    comp_name = 'HSLA (Halosteric)';
    if exist('sig_HSLA', 'var'), sig_data = sig_HSLA; else, sig_data = zeros(size(trend_data)); end
    fprintf('>> ✅ 成功提取并锁定盐收缩 HSLA 的趋势数据。\n');

else
    error('❌ 无效的目标变量: %s。请使用 ''SSLA'', ''TSLA'' 或 ''HSLA''', target_comp);
end

% ！！！智能单位换算 (m -> mm) ！！！
if max(abs(trend_data(:))) < 0.2
    fprintf('>> 🔄 检测到单位可能为 [米(m)]，自动乘以 1000 转换为 [毫米(mm)]...\n');
    trend_data = trend_data * 1000;
end

% 经纬度处理与对齐
if ~exist('Lon', 'var') && exist('target_lon', 'var'), Lon = target_lon; end
if ~exist('Lat', 'var') && exist('target_lat', 'var'), Lat = target_lat; end
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% 第一步：先转置使得矩阵符合地图绘制要求 [Lat x Lon]
if size(trend_data, 1) ~= length(Lat)
    temp_t = nan(length(Lat), length(Lon), size(trend_data, 3));
    temp_s = nan(length(Lat), length(Lon), size(trend_data, 3));
    for i = 1:size(trend_data, 3)
        temp_t(:,:,i) = trend_data(:,:,i)';
        temp_s(:,:,i) = sig_data(:,:,i)';
    end
    trend_data = temp_t;
    sig_data = temp_s;
end

% 第二步：再进行经度 [0, 360] -> [-180, 180] 转换与切片拼接
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    for i = 1:size(trend_data, 3)
        trend_data(:,:,i) = trend_data(:, sort_idx, i);
        sig_data(:,:,i) = sig_data(:, sort_idx, i);
    end
end
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

%% 3. 绘制空间分布图 (顶刊 RdBu 色棒)
fprintf('>> 绘制趋势空间分布图...\n');
max_terms = size(trend_data, 3);
if max_terms > 9, max_terms = 9; end % 扩容支持 9 阶

figure('Position', [100, 50, 1200, max(400, ceil(max_terms/2)*350)], 'Color', 'w', 'Name', 'Spatial Maps');
t = tiledlayout(ceil(max_terms/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    data = trend_data(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    % 显著性打点：在不显著 (== 0) 的区域打深灰点掩盖
    sig_mask = sig_data(:,:,term_idx); skip = 4; 
    [row, ~] = find(sig_mask(1:skip:end, 1:skip:end) == 0);
    if ~isempty(row)
        sig_lon = Lon_Grid(1:skip:end, 1:skip:end); sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==0), sig_lat(sig_mask(1:skip:end, 1:skip:end)==0), ...
               'o', 'markerfacecolor', [.65 .65 .65], 'markeredgecolor', 'none', 'markersize', 0.8);
    end
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap);
    
    if strcmp(state, 'Average')
        if term_idx == 1, pct = 99.5; elseif term_idx <= 3, pct = 99.0; else, pct = 98.8; end
    else
        if term_idx <= 2, pct = 99.5; elseif term_idx <= 4, pct = 99.0; else, pct = 98.5; end
    end
    
    raw_max = prctile(abs(data(:)), pct);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([-final_limit, final_limit]); 
    
    % 🌟 智能子图标题
    if max_terms == 9 && term_idx == 1
        title('Exact SSLA (Direct Difference)', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1]);
    elseif max_terms == 9
        title(sprintf('Taylor Approx Order %d', term_idx-1), 'FontSize', 12, 'FontWeight', 'bold');
    else
        title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end
sgtitle(sprintf('%s (%s) Trend Distribution [%s]', dataset_name, state, strtok(comp_name, ' ')), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 4. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(trend_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Trend = nan(num_valid, max_terms);
for n = 1:max_terms
    data_slice = trend_data(:,:,n);
    vec_t = nan(num_valid, 1);
    for k = 1:num_valid
        vec_t(k) = data_slice(rows(k), cols(k));
    end
    Data_Trend(:,n) = vec_t;
end

%% 5. 绘制箱线图 (🚨 修复 X轴着色丢失与 Y轴挤压问题)
fprintf('>> 绘制趋势箱线图...\n');
figure('Position', [50, 50, 1300, 550], 'Color', 'w', 'Name', 'Trend Boxplot');

C1 = [0.902, 0.294, 0.208]; C2 = [0.302, 0.455, 0.690]; C3 = [0.000, 0.627, 0.408]; C4 = [0.5 0 0.5]; C0 = [0.3 0.3 0.3];
BoxBgColor = [0.96 0.96 0.96]; BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.25; ScatterSize = 10;
X_Limit = [0.4, max_terms + 1.5]; 
% 预留更大的右侧空间防止轴被吃掉
MainPos = [0.06, 0.12, 0.72, 0.73]; 

% =======================================================
% 🎯 动态自适应分组引擎
% =======================================================
if max_terms == 9 % SSLA 合成版 (1层黑 + 8层彩)
    if strcmp(state, 'Average')
        G_idx = { [1, 2], [3, 4], [5, 6], [7, 8, 9] };
        Pos_G = { [1, 2.2], [3.4, 4.4], [5.6, 6.6], [7.8, 8.8, 9.8] };
    else
        G_idx = { [1, 2, 3, 4], [5, 6], [7, 8, 9] };
        Pos_G = { [1, 2.2, 3.2, 4.2], [5.4, 6.4], [7.6, 8.6, 9.6] };
    end
    c_list = {C1, C2, C3, C4};
    x_tick_labels = {};
    for i=1:length(G_idx)
        for j=1:length(G_idx{i})
            k = G_idx{i}(j);
            if k == 1, c_str = sprintf('\\color[rgb]{%f,%f,%f}', C0); lbl = 'Exact';
            else, c_str = sprintf('\\color[rgb]{%f,%f,%f}', c_list{i}); lbl = sprintf('O%d', k-1); end
            x_tick_labels{end+1} = [c_str lbl];
        end
    end

elseif max_terms == 8 % TSLA / HSLA 标准版 (纯 8 层)
    if strcmp(state, 'Average')
        G_idx = { [1], [2, 3], [4, 5], [6, 7, 8] };
        Pos_G = { [1], [2.4, 3.4], [4.8, 5.8], [7.2, 8.2, 9.2] };
    else
        G_idx = { [1, 2, 3], [4, 5, 6], [7, 8] };
        Pos_G = { [1, 2, 3], [4.4, 5.4, 6.4], [7.8, 8.8] };
    end
    c_list = {C1, C2, C3, C4};
    x_tick_labels = {};
    for i=1:length(G_idx)
        for j=1:length(G_idx{i})
            k = G_idx{i}(j);
            c_str = sprintf('\\color[rgb]{%f,%f,%f}', c_list{i}); lbl = sprintf('T%d', k);
            x_tick_labels{end+1} = [c_str lbl];
        end
    end

else % 兜底 (比如只有 1 层)
    G_idx = { 1:max_terms };
    Pos_G = { 1:max_terms };
    c_list = {C1};
    x_tick_labels = {sprintf('\\color[rgb]{%f,%f,%f}Total Trend', C1)};
end

% 开始画坐标轴层
num_layers = length(G_idx);
Axes_List = gobjects(1, num_layers);
Scale_Factors = ones(1, num_layers);
Lim_Min = zeros(1, num_layers); Lim_Max = zeros(1, num_layers);

for i = 1:num_layers
    curr_data = Data_Trend(:, G_idx{i});
    max_v = max(curr_data(:)); min_v = min(curr_data(:));
    range_v = max_v - min_v; if range_v < 1e-12, range_v = 2; min_v = -1; max_v = 1; end
    Lim_Min(i) = min_v - range_v * 0.1; Lim_Max(i) = max_v + range_v * 0.1;
    
    if i == 1
        ax = axes('Position', MainPos);
    else
        % 🚨 核心修复：大幅缩小偏移比例，防止右侧坐标轴被挤出画布
        Offset_Ratio = 0.05 * (i-1);
        Pos_curr = MainPos; Pos_curr(3) = MainPos(3) * (1 + Offset_Ratio);
        ax = axes('Position', Pos_curr);
        Scale_Factors(i) = Pos_curr(3) / MainPos(3);
    end
    hold(ax, 'on'); Axes_List(i) = ax;
    
    for j = 1:length(G_idx{i})
        k = G_idx{i}(j);
        if max_terms == 9 && k == 1, curr_c = C0; else, curr_c = c_list{i}; end
        draw_jitter_unit(ax, Pos_G{i}(j), Data_Trend(:,k), BoxWidth, ScatterWidth, ScatterSize, curr_c, BoxBgColor, ScatterAlpha);
    end
    
    ylim(ax, [Lim_Min(i), Lim_Max(i)]);
    if i == 1
        
        ylabel(ax, 'Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
        % 第一层坐标轴颜色如果只有一层设为主色，否则为黑色
        if num_layers == 1, y_col = C1; else, y_col = 'k'; end
        set(ax, 'XLim', X_Limit, 'YColor', y_col, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    else
        Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factors(i)];
        set(ax, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', c_list{i}, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    end
end

% 🚨 核心修复：重新挂载强力 X 轴刻度及彩色标签
All_Ticks = [];
for i = 1:num_layers
    All_Ticks = [All_Ticks, Pos_G{i}];
    if i < num_layers
        xline(Axes_List(1), mean([Pos_G{i}(end), Pos_G{i+1}(1)]), ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    end
end
% 强制激活 tex 解析器渲染颜色
set(Axes_List(1), 'XTick', All_Ticks, 'XTickLabel', x_tick_labels, 'FontWeight', 'bold', 'TickLabelInterpreter', 'tex');

% ==========================================================
% 统一通用图例与标题 
% ==========================================================
xlabel(Axes_List(1), sprintf('%s Components', strtok(comp_name, ' ')), 'FontSize', 16, 'FontWeight', 'bold');
ax_title = axes('Position', [0, 0.92, 1, 0.05], 'Visible', 'off');
text(ax_title, 0.5, 0.5, sprintf('%s (%s) %s Trend Boxplot', dataset_name, state, strtok(comp_name, ' ')), 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');

% 🚨 调整横向图例的起始位置，避免文字超出画布
ax_lgd = axes('Position', [MainPos(1)-0.02, 0.88, MainPos(3)+0.1, 0.04], 'XLim', [0, 100], 'YLim', [0, 1], 'Visible', 'off'); hold(ax_lgd, 'on');
patch(ax_lgd, [3, 5, 5, 3], [0.3, 0.3, 0.7, 0.7], BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
text(ax_lgd, 6, 0.5, 'IQR Box', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [19, 22], [0.5, 0.5], 'r-', 'LineWidth', 2);
text(ax_lgd, 23, 0.5, 'Median line', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [37], [0.5], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
text(ax_lgd, 39, 0.5, 'Mean value', 'FontSize', 13, 'FontName', 'Times New Roman');

plot(ax_lgd, [55, 58], [0.5, 0.5], 'k-', 'LineWidth', 1.5);
text(ax_lgd, 59, 0.5, '2.5% \sim 97.5%', 'Interpreter', 'tex', 'FontSize', 13, 'FontName', 'Times New Roman');

scatter(ax_lgd, 78, 0.5, 50, C1, 'filled', 'MarkerFaceAlpha', 0.8);
text(ax_lgd, 80, 0.5, sprintf('Outliers (%s)', dataset_name), 'FontSize', 13, 'FontName', 'Times New Roman');

if num_layers > 1
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on'; hPan = pan(gcf); hPan.Motion = 'both'; 
    sync_cb = @(s,e) SyncDynAxes(e, Axes_List, Lim_Min, Lim_Max, Scale_Factors);
    hZoom.ActionPostCallback = sync_cb; hPan.ActionPostCallback = sync_cb;
end

fprintf('>> 🎉 完美出图！\n');

% 检查批处理环境并保存
save_dir = 'D:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i); if isvalid(fig)
        fig_name = strrep(strrep(get(fig, 'Name'), ' ', '_'), '(', ''); fig_name = strrep(fig_name, ')', '');
        filename = fullfile(save_dir, sprintf('Trend_%s_%s_%s_%s.jpg', dataset_name, state, strtok(comp_name, ' '), fig_name));
        try exportgraphics(fig, filename, 'Resolution', 300); catch, saveas(fig, filename, 'jpg'); end
    end
end
end

%% ============================================================
%% 内部函数
%% ============================================================
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); if isempty(data), return; end
    q1 = prctile(data, 25); q3 = prctile(data, 75); med_val = median(data); mean_val = mean(data); low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    
    idx_out = data < low_w | data > up_w; outliers = data(idx_out);
    if ~isempty(outliers)
        if length(outliers) > 3000, outliers = outliers(randperm(length(outliers), 3000)); end
        x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
        try scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
        catch, scatter(ax, x_jit, outliers, sz_scatter, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s); end
    end
    
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2); plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    cap_w = w_box * 0.3; plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2); plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    x_L = center - w_box/2; x_R = center + w_box/2; patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2); plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

function SyncDynAxes(event_obj, Axes_List, Lim_Min, Lim_Max, Scale_Factors)
    try
        curr_ax = event_obj.Axes; if isempty(curr_ax), curr_ax = Axes_List(1); end
        idx = find(Axes_List == curr_ax, 1); if isempty(idx), idx = 1; end
        
        y_curr = get(curr_ax, 'YLim');
        curr_range = Lim_Max(idx) - Lim_Min(idx);
        ratio_min = (y_curr(1) - Lim_Min(idx)) / curr_range;
        ratio_max = (y_curr(2) - Lim_Min(idx)) / curr_range;
        
        curr_x_lim = get(curr_ax, 'XLim');
        x1_base = curr_x_lim(1);
        x1_range = (curr_x_lim(2) - curr_x_lim(1)) / Scale_Factors(idx);
        
        for i = 1:length(Axes_List)
            i_range = Lim_Max(i) - Lim_Min(i);
            set(Axes_List(i), 'YLim', [Lim_Min(i) + ratio_min * i_range, Lim_Min(i) + ratio_max * i_range]);
            if i == 1
                set(Axes_List(i), 'XLim', [x1_base, x1_base + x1_range]);
            else
                set(Axes_List(i), 'XLim', [x1_base, x1_base + x1_range * Scale_Factors(i)]);
            end
        end
    catch; end
end