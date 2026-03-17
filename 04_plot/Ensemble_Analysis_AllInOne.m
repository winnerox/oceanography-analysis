function Ensemble_Analysis_AllInOne(state)
% 功能: 多模式集合平均与置信区间分析，支持空间插值、T检验及多源溯源箱线图
% 特性: 自动量级聚类分层、修复X轴漂移对齐Bug、科学防过曝、高透星云散点、清晰显著性打点

clc; close all;
fprintf('\n======================================================\n');
fprintf('🚀 启动多源集合分析引擎 (完美对齐 & 清晰打点版) [%s]\n', state);
fprintf('======================================================\n');

%% 0. 核心控制开关 (🚨 绝对锁定目标变量)
comp_type = 'TSLA';  % 可选: 'TSLA' (热膨胀) 或 'HSLA' (盐收缩) 或 'SSLA' (总比容)
fprintf('🎯 当前锁定分析目标: 【%s】\n', comp_type);

%% 1. 数据加载与空间插值对齐
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
amp_dir = 'D:\work\Task_Convergence\Amplitude_Results';
datasets = {'EN4', 'IAP', 'Ishii'};

Trend_All = []; Amp_All = [];

% === 设定标准大洋网格 ===
target_lon = 0.5:1:359.5; target_lat = -89.5:1:89.5;
[Target_Lon_Grid, Target_Lat_Grid] = meshgrid(target_lon, target_lat);

for d = 1:3
    trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', datasets{d}, state));
    amp_file = fullfile(amp_dir, sprintf('%s_%s_Amplitude.mat', datasets{d}, state));
    if ~exist(trend_file, 'file'), error('找不到趋势文件: %s', trend_file); end
    if ~exist(amp_file, 'file'), error('找不到振幅文件: %s', amp_file); end
    
    t_data = load(trend_file); a_data = load(amp_file);
    
    var_trend = sprintf('trend_%s', comp_type); var_amp = sprintf('amp_%s', comp_type);
    if ~isfield(t_data, var_trend) || isempty(t_data.(var_trend))
        error('❌ 在 %s 中找不到指定的变量: %s', datasets{d}, var_trend);
    end
    
    T_curr = t_data.(var_trend); A_curr = abs(a_data.(var_amp)); 
    curr_lon = double(t_data.Lon(:))'; curr_lat = double(t_data.Lat(:))';
    
    n_lon = length(curr_lon); n_lat = length(curr_lat);
    if size(T_curr, 1) == n_lon && size(T_curr, 2) == n_lat
        T_curr = permute(T_curr, [2, 1, 3]); A_curr = permute(A_curr, [2, 1, 3]);
    end
    
    % 二维空间插值对齐 (带周期性边界)
    need_interp = false;
    if n_lat ~= 180 || n_lon ~= 360, need_interp = true;
    elseif abs(curr_lon(1) - target_lon(1)) > 0.1 || abs(curr_lon(end) - target_lon(end)) > 0.1, need_interp = true; end
    
    if need_interp
        if min(curr_lon) < 0
            curr_lon(curr_lon < 0) = curr_lon(curr_lon < 0) + 360;
            [curr_lon, sort_idx] = sort(curr_lon);
            T_curr = T_curr(:, sort_idx, :); A_curr = A_curr(:, sort_idx, :);
        end
        curr_lon = [curr_lon(end)-360, curr_lon, curr_lon(1)+360];
        T_curr = cat(2, T_curr(:,end,:), T_curr, T_curr(:,1,:)); A_curr = cat(2, A_curr(:,end,:), A_curr, A_curr(:,1,:));
        [Curr_Lon_Grid, Curr_Lat_Grid] = meshgrid(curr_lon, curr_lat);
        T_interp = nan(180, 360, size(T_curr,3)); A_interp = nan(180, 360, size(A_curr,3));
        for k = 1:size(T_curr, 3)
            T_interp(:,:,k) = interp2(Curr_Lon_Grid, Curr_Lat_Grid, T_curr(:,:,k), Target_Lon_Grid, Target_Lat_Grid, 'linear');
            A_interp(:,:,k) = interp2(Curr_Lon_Grid, Curr_Lat_Grid, A_curr(:,:,k), Target_Lon_Grid, Target_Lat_Grid, 'linear');
        end
        T_curr = T_interp; A_curr = A_interp;
    end
    Trend_All = cat(4, Trend_All, T_curr); Amp_All = cat(4, Amp_All, A_curr);
end

Lon = target_lon; Lat = target_lat;
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    Trend_All = Trend_All(:, sort_idx, :, :); Amp_All = Amp_All(:, sort_idx, :, :);
end
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);
max_terms = size(Trend_All, 3); if max_terms > 8, max_terms = 8; end

%% 2. 计算集合平均与 90% 显著性检验
fprintf('>> 计算集合平均与 90%% 置信区间...\n');
Trend_Ens = mean(Trend_All, 4, 'omitnan'); Amp_Ens = mean(Amp_All, 4, 'omitnan');
n_sample = 3; df = n_sample - 1; alpha = 0.10;             
std_field = std(Trend_All, 0, 4, 'omitnan');         
standard_error = std_field ./ sqrt(n_sample);        
t_stat = Trend_Ens ./ standard_error;                
P_value = 2 * (1 - tcdf(abs(t_stat), df)); 
Sig_Ens = P_value < alpha; % P < 0.1 判定为显著 (1)

%% 3. 准备箱线图原始数据
Box_Trend_Raw = cell(3, max_terms); Box_Amp_Raw = cell(3, max_terms); 
for n = 1:max_terms
    for d = 1:3
        raw_T = Trend_All(:,:,n,d); Box_Trend_Raw{d,n} = raw_T(~isnan(raw_T));
        raw_A = Amp_All(:,:,n,d);   Box_Amp_Raw{d,n} = raw_A(~isnan(raw_A));
    end
end

%% 4. 色盘与颜色定义
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
Cmap_Trend = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];
Amp_Points = [1.00 1.00 1.00; 1.00 0.90 0.50; 0.95 0.50 0.15; 0.80 0.10 0.15; 0.40 0.00 0.00]; Amp_Locs = [0, 0.2, 0.5, 0.8, 1]; 
Cmap_Amp = [interp1(Amp_Locs, Amp_Points(:,1), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,2), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,3), xx, 'pchip')'];

% EN4(红), IAP(蓝), Ishii(绿)
colors_raw = {[0.902, 0.294, 0.208], [0.302, 0.455, 0.690], [0.000, 0.627, 0.408]}; 
BoxBgColor = [0.96 0.96 0.96]; 

%% 5. 绘制图表
fprintf('>> 开始绘制顶刊图表...\n');
Draw_Spatial_Map(Trend_Ens, Sig_Ens, max_terms, Lon_Grid, Lat_Grid, Cmap_Trend, sprintf('Ensemble Mean %s Trend (%s)', comp_type, state), 'mm/yr', true, state);
Draw_Spatial_Map(Amp_Ens, [], max_terms, Lon_Grid, Lat_Grid, Cmap_Amp, sprintf('Ensemble Mean %s Amplitude (%s)', comp_type, state), 'mm', false, state);

if max_terms > 1
    Draw_Grouped_MultiAxis_Boxplot(max_terms, Box_Trend_Raw, state, colors_raw, BoxBgColor, sprintf('Ensemble %s Trend Comparison (%s)', comp_type, state), 'Trend (mm/yr)', true);
    Draw_Grouped_MultiAxis_Boxplot(max_terms, Box_Amp_Raw, state, colors_raw, BoxBgColor, sprintf('Ensemble %s Amplitude Comparison (%s)', comp_type, state), 'Amplitude (mm)', false);
end

%% 6. 保存
save_dir = 'D:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i); if isvalid(fig)
        fig_name = strrep(strrep(strrep(get(fig, 'Name'), ' ', '_'), '(', ''), ')', '');
        filename = fullfile(save_dir, sprintf('Ensemble_%s_%s.jpg', state, fig_name));
        try exportgraphics(fig, filename, 'Resolution', 300); catch, saveas(fig, filename, 'jpg'); end
        fprintf('   [成功] 已导出: %s\n', filename);
    end
end
fprintf('🎉 恭喜！终极分析完毕！\n');
end

%% =========================================================================
%% 内部函数 1: 顶刊地图绘制 (🚨 强化显著性打点视觉效果)
function Draw_Spatial_Map(Data_Map, Sig_Map, max_terms, Lon_Grid, Lat_Grid, cmap, fig_title, unit_str, is_divergent, state)
    layout_rows = ceil(max_terms / 2);
    figure('Position', [100, 50, 1200, max(400, layout_rows*350)], 'Color', 'w', 'Name', fig_title);
    t = tiledlayout(layout_rows, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for i = 1:max_terms
        nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        data = Data_Map(:,:,i); m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
        
        % =========================================================
        % 🌟 核心优化：高可见度显著性打点 (Stippling)
        % =========================================================
        if ~isempty(Sig_Map)
            sig_mask = Sig_Map(:,:,i); 
            skip = 3; % 缩小间隔，让点更密集 (原来是6)
            
            % 提取不显著区域 (0 代表不显著)
            [row, ~] = find(sig_mask(1:skip:end, 1:skip:end) == 0); 
            if ~isempty(row)
                sig_lon = Lon_Grid(1:skip:end, 1:skip:end); 
                sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
                
                % 采用浅灰色和小尺寸，以免喧宾夺主
                m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==0), ...
                       sig_lat(sig_mask(1:skip:end, 1:skip:end)==0), ...
                       'o', 'markerfacecolor', [.65 .65 .65], 'markeredgecolor', 'none', 'markersize', 0.8);
            end
        end
        
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); colormap(gca, cmap);
        
        if strcmp(state, 'Average')
            if i == 1, pct = 99.5; elseif i <= 3, pct = 99.0; else, pct = 98.8; end
        else
            if i <= 2, pct = 99.5; elseif i <= 4, pct = 99.0; else, pct = 98.5; end
        end
        
        if is_divergent
            raw_max = prctile(abs(data(:)), pct); 
        else
            raw_max = prctile(data(:), pct); 
        end
        if isempty(raw_max) || isnan(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
        
        exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
        if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
        final_limit = nice_base * 10^exponent;
        
        if is_divergent, caxis([-final_limit, final_limit]); else, caxis([0, final_limit]); end
        title(sprintf('Term %d', i), 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
        title(cb, unit_str, 'FontSize', 8, 'FontWeight', 'bold'); 
        if is_divergent, cb.Ticks = linspace(-final_limit, final_limit, 5); else, cb.Ticks = linspace(0, final_limit, 5); end
    end
    sgtitle(fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
end

%% 内部函数 2: 🚀 终极动态分层箱线图 (自动量级聚类 + 完美坐标轴对齐)
function Draw_Grouped_MultiAxis_Boxplot(max_terms, Raw_Cell, state, colors_raw, BoxBgColor, fig_title, y_label, show_zero_band)
    figure('Position', [50, 50, 1300, 580], 'Color', 'w', 'Name', fig_title);
    
    AxC = {[0.1 0.1 0.1], [0.5 0 0.5], [0.8 0.4 0], [0 0.5 0.6]};
    BoxWidth = 0.55; ScatterWidth = 0.4; ScatterAlpha = 0.20; ScatterSize = 10;
    MainPos = [0.07, 0.12, 0.76, 0.72]; 
    
    mags = zeros(1, max_terms);
    for k = 1:max_terms
        temp_data = [];
        for d = 1:3
            vec = Raw_Cell{d,k}(:); temp_data = [temp_data; vec(~isnan(vec))];
        end
        if isempty(temp_data), mags(k) = -12; continue; end
        val = prctile(abs(temp_data), 98); if val < 1e-12, val = 1e-12; end
        mags(k) = log10(val); 
    end
    
    diffs = mags(1:end-1) - mags(2:end); 
    [sort_diffs, sort_idx] = sort(diffs, 'descend');
    
    valid_splits = sort_idx(sort_diffs > 0.45); 
    num_splits = min(3, length(valid_splits)); 
    split_points = sort(valid_splits(1:num_splits));
    
    G_idx = cell(1, length(split_points) + 1); Pos_G = cell(1, length(split_points) + 1);
    start_idx = 1; current_x = 1;
    for i = 1:length(split_points)
        G_idx{i} = start_idx : split_points(i);
        Pos_G{i} = current_x : 1 : (current_x + length(G_idx{i}) - 1);
        current_x = Pos_G{i}(end) + 1.4; 
        start_idx = split_points(i) + 1;
    end
    G_idx{end} = start_idx : max_terms;
    Pos_G{end} = current_x : 1 : (current_x + length(G_idx{end}) - 1);
    
    num_layers = length(G_idx);
    X_Limit = [Pos_G{1}(1) - 0.6, Pos_G{end}(end) + 1.0];
    
    Lim_Min = zeros(1, num_layers); Lim_Max = zeros(1, num_layers);
    for i = 1:num_layers
        d_curr = [];
        for k = G_idx{i}
            for d = 1:3
                tmp = Raw_Cell{d,k}(:); if ~show_zero_band, tmp = abs(tmp); end
                d_curr = [d_curr; tmp(~isnan(tmp))];
            end
        end
        if show_zero_band
            min_v = min(d_curr); max_v = max(d_curr); 
            if isnan(min_v)||isnan(max_v)||(max_v-min_v)<1e-12, min_v=-1; max_v=1; end
            range = max_v - min_v; Lim_Min(i) = min_v - range*0.1; Lim_Max(i) = max_v + range*0.1;
        else
            max_v = max(d_curr); if isnan(max_v)||max_v==0, max_v=1; end
            Lim_Min(i) = 0; Lim_Max(i) = max_v * 1.1;
        end
    end
    Axes_List = gobjects(1, num_layers);
    Scale_Factors = ones(1, num_layers);
    
    for i = 1:num_layers
        if i == 1
            ax = axes('Position', MainPos);
        else
            Offset_Ratio = 0.055 * (i-1);
            Pos_curr = MainPos; Pos_curr(3) = MainPos(3) * (1 + Offset_Ratio);
            ax = axes('Position', Pos_curr);
            Scale_Factors(i) = Pos_curr(3) / MainPos(3);
        end
        hold(ax, 'on'); Axes_List(i) = ax;
        
        for j = 1:length(G_idx{i})
            k = G_idx{i}(j); data_c = cell(1,3); 
            for d=1:3, data_c{d}=Raw_Cell{d,k}; if ~show_zero_band, data_c{d}=abs(data_c{d}); end; end
            draw_colored_jitter_unit(ax, Pos_G{i}(j), data_c, BoxWidth, ScatterWidth, ScatterSize, colors_raw, BoxBgColor, ScatterAlpha);
        end
        
        ylim(ax, [Lim_Min(i), Lim_Max(i)]);
        if i == 1
            if show_zero_band, yline(ax, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.2); 
            else, yline(ax, 0, '-', 'Color', 'k', 'LineWidth', 1.5); end
            ylabel(ax, y_label, 'FontWeight', 'bold', 'FontSize', 16);
            set(ax, 'XLim', X_Limit, 'YColor', AxC{i}, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        else
            Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factors(i)];
            set(ax, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', AxC{i}, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        end
    end
    All_Ticks = []; All_Labels = {};
    for i = 1:num_layers
        All_Ticks = [All_Ticks, Pos_G{i}];
        for k = G_idx{i}, All_Labels{end+1} = sprintf('\\color[rgb]{%f,%f,%f}T%d', AxC{i}, k); end
        if i < num_layers
            xline(Axes_List(1), mean([Pos_G{i}(end), Pos_G{i+1}(1)]), ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
        end
    end
    
    set(Axes_List(1), 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold', 'TickLabelInterpreter', 'tex');
    xlabel(Axes_List(1), 'Taylor Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    
    ax_title = axes('Position', [0, 0.93, 1, 0.05], 'Visible', 'off');
    text(ax_title, 0.5, 0.5, fig_title, 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');
    
    ax_lgd = axes('Position', [MainPos(1), 0.88, MainPos(3), 0.04], 'XLim', [0, 100], 'YLim', [0, 1], 'Visible', 'off'); hold(ax_lgd, 'on');
    plot(ax_lgd, [5, 10], [0.5, 0.5], 'r-', 'LineWidth', 2); text(ax_lgd, 11, 0.5, 'Median line', 'FontSize', 13, 'FontName', 'Times New Roman');
    plot(ax_lgd, [23, 28], [0.5, 0.5], 'k-', 'LineWidth', 1.5); text(ax_lgd, 29, 0.5, '2.5% \sim 97.5%', 'Interpreter', 'tex', 'FontSize', 13, 'FontName', 'Times New Roman');
    
    scatter(ax_lgd, 47, 0.5, 60, colors_raw{1}, 'filled', 'MarkerFaceAlpha', 0.9); text(ax_lgd, 49, 0.5, 'EN4 Outliers', 'FontSize', 13, 'FontName', 'Times New Roman');
    scatter(ax_lgd, 65, 0.5, 60, colors_raw{2}, 'filled', 'MarkerFaceAlpha', 0.9); text(ax_lgd, 67, 0.5, 'IAP Outliers', 'FontSize', 13, 'FontName', 'Times New Roman');
    scatter(ax_lgd, 83, 0.5, 60, colors_raw{3}, 'filled', 'MarkerFaceAlpha', 0.9); text(ax_lgd, 85, 0.5, 'Ishii Outliers', 'FontSize', 13, 'FontName', 'Times New Roman');
    
    if num_layers > 1
        hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on'; hPan = pan(gcf); hPan.Motion = 'both'; 
        sync_cb = @(s,e) SyncDynAxes(e, Axes_List, Lim_Min, Lim_Max, Scale_Factors, show_zero_band);
        hZoom.ActionPostCallback = sync_cb; hPan.ActionPostCallback = sync_cb;
    end
end

%% 内部核心引擎: 绘制统一箱体 + 彩色溯源实心球体
function draw_colored_jitter_unit(ax, center, data_cell, w_box, w_scatter, sz_scatter, colors_raw, color_box, alpha_s)
    data_all = []; for d = 1:3, data_all = [data_all; data_cell{d}(~isnan(data_cell{d}))]; end
    if isempty(data_all), return; end
    
    q1 = prctile(data_all, 25); q3 = prctile(data_all, 75);
    med_val = median(data_all); mean_val = mean(data_all);
    low_w = prctile(data_all, 2.5); up_w  = prctile(data_all, 97.5);
    
    for d = 1:3
        data_d = data_cell{d}(~isnan(data_cell{d})); outliers = data_d(data_d < low_w | data_d > up_w);
        if ~isempty(outliers)
            if length(outliers) > 3000, outliers = outliers(randperm(length(outliers), 3000)); end
            x_jit = center + (rand(size(outliers)) - 0.5) * w_scatter;
            color_dots = colors_raw{d}; 
            try scatter(ax, x_jit, outliers, sz_scatter * 1.5, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s);
            catch, scatter(ax, x_jit, outliers, sz_scatter * 1.5, color_dots, 'filled', 'MarkerFaceAlpha', alpha_s); end
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
end

%% 内部辅助函数: 全局任意自适应层级的联动同步引擎
function SyncDynAxes(event_obj, Axes_List, Lim_Min, Lim_Max, Scale_Factors, has_zero)
    try
        curr_ax = event_obj.Axes; if isempty(curr_ax), curr_ax = Axes_List(1); end
        idx = find(Axes_List == curr_ax, 1); if isempty(idx), idx = 1; end
        
        y_curr = get(curr_ax, 'YLim');
        if has_zero
            curr_range = Lim_Max(idx) - Lim_Min(idx);
            ratio_min = (y_curr(1) - Lim_Min(idx)) / curr_range;
            ratio_max = (y_curr(2) - Lim_Min(idx)) / curr_range;
        else
            ratio_min = 0; ratio_max = y_curr(2) / Lim_Max(idx);
        end
        
        curr_x_lim = get(curr_ax, 'XLim');
        x1_base = curr_x_lim(1);
        x1_range = (curr_x_lim(2) - curr_x_lim(1)) / Scale_Factors(idx);
        
        for i = 1:length(Axes_List)
            if has_zero
                i_range = Lim_Max(i) - Lim_Min(i);
                set(Axes_List(i), 'YLim', [Lim_Min(i) + ratio_min * i_range, Lim_Min(i) + ratio_max * i_range]);
            else
                set(Axes_List(i), 'YLim', [0, Lim_Max(i) * ratio_max]);
            end
            if i == 1
                set(Axes_List(i), 'XLim', [x1_base, x1_base + x1_range]);
            else
                set(Axes_List(i), 'XLim', [x1_base, x1_base + x1_range * Scale_Factors(i)]);
            end
        end
    catch; end
end