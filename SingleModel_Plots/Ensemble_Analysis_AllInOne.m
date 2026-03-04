function Ensemble_Analysis_AllInOne(state)
% 功能: 一键完成三套数据的集合平均、极速统计计算，并绘制4张顶刊质量图表
% 特色: 保留稳健的数据维度读取，采用“共用X轴 + 统一灰箱体 + 溯源彩色散点 + 多Y轴”的终极设计
% 参数: state = 'StdRef' 或 'Average'

clc; close all;
fprintf('\n======================================================\n');
fprintf('🚀 启动集合分析引擎 (统一X轴溯源散点版) [%s]\n', state);
fprintf('======================================================\n');

%% 1. 数据加载与集合平均计算 (完美保留你的读取对齐逻辑)
datasets = {'EN4', 'IAP', 'Ishii'};
data_dirs = {
    'D:\work\EN4_TSLA_Terms';
    'D:\work\IAP_mat_data';
    'D:\work\Ishii_mat_data'
};
TSLA_All = [];

dataset_sizes = [];
dataset_data = cell(1, 3);
dataset_coords = cell(1, 3);

for i = 1:3
    DataDir = data_dirs{i};
    file_name = fullfile(DataDir, sprintf('%s_TSLA_Terms_1to8_%s.mat', datasets{i}, state));
    if ~exist(file_name, 'file')
        file_name = sprintf('%s_TSLA_Terms_1to8_%s.mat', datasets{i}, state); 
    end
    fprintf('>> 加载 %s 数据...\n', datasets{i});
    data = load(file_name);
    
    if isfield(data, 'TSLA_AllOrders'), tsla_raw = data.TSLA_AllOrders;
    elseif isfield(data, 'TSLA_Result'), tsla_raw = data.TSLA_Result; end
    
    if isfield(data, 'lon'), lon = data.lon; elseif isfield(data, 'Lon'), lon = data.Lon; end
    if isfield(data, 'lat'), lat = data.lat; elseif isfield(data, 'Lat'), lat = data.Lat; end
    if isfield(data, 'time_axis'), t_vec = data.time_axis(:); elseif isfield(data, 'Time_Axis'), t_vec = data.Time_Axis(:); end
    
    dataset_data{i} = tsla_raw;
    dataset_coords{i} = struct('lon', lon, 'lat', lat, 'time', t_vec);
    dataset_sizes(i,:) = size(tsla_raw);
end

[~, min_idx] = min(dataset_sizes(:,2)); 
base_size = dataset_sizes(min_idx,:);
base_coords = dataset_coords{min_idx};
lon = base_coords.lon; lat = base_coords.lat; t_vec = base_coords.time;

TSLA_All = zeros([base_size, 3]);

for i = 1:3
    tsla_raw = dataset_data{i};
    current_size = size(tsla_raw);
    if ~isequal(current_size, base_size)
        fprintf('   调整 %s 数据维度以匹配基准...\n', datasets{i});
        if current_size(2) ~= base_size(2)
            if current_size(2) > base_size(2)
                tsla_raw = tsla_raw(:, 1:base_size(2), :, :);
            else
                padding = base_size(2) - current_size(2);
                tsla_raw = cat(2, tsla_raw, nan(size(tsla_raw,1), padding, size(tsla_raw,3), size(tsla_raw,4)));
            end
        end
    end
    TSLA_All(:,:,:,:,i) = tsla_raw;
end

fprintf('>> 计算网格级集合平均 (Ensemble Mean)...\n');
TSLA_Ens = mean(TSLA_All, 5, 'omitnan'); 

if max(lon) > 180
    lon(lon > 180) = lon(lon > 180) - 360;
    [lon, sort_idx] = sort(lon);
    TSLA_Ens = TSLA_Ens(sort_idx, :, :, :);
    TSLA_All = TSLA_All(sort_idx, :, :, :, :);
end
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);
nlon = length(lon); nlat = length(lat); ntime = length(t_vec);

%% 2. 向量化极速计算 (不仅算平均，还算三家独立数据以提取散点)
fprintf('>> 向量化极速计算趋势、显著性与振幅 (耗时 < 1秒)...\n');
Map_Trend = nan(nlat, nlon, 8); Map_Sig = nan(nlat, nlon, 8); Map_Amp = nan(nlat, nlon, 8);
Box_Trend_Ens = cell(1, 8); Box_Amp_Ens = cell(1, 8); 
Box_Trend_Raw = cell(3, 8); Box_Amp_Raw = cell(3, 8); 

X_mat = [t_vec - mean(t_vec), ones(ntime, 1)];
invXTX = (X_mat' * X_mat) \ eye(2);

for n = 1:8
    % --- 1. 处理集合平均 (生成地图和中央灰箱子) ---
    slice_ens = squeeze(TSLA_Ens(:,:,:,n)); 
    Y_ens = reshape(slice_ens, nlon*nlat, ntime)'; 
    valid_idx = ~isnan(Y_ens(1,:));
    Y_valid = Y_ens(:, valid_idx);
    
    B = X_mat \ Y_valid; 
    trend_valid = B(1,:) * 1000; 
    Resid = Y_valid - X_mat * B;
    amp_valid = std(Resid, 0, 1) * 1000; 
    
    sigma_sq = sum(Resid.^2, 1) / (ntime - 2);
    se_slope = sqrt(sigma_sq * invXTX(1,1));
    t_stat = B(1,:) ./ se_slope;
    p_val = 2 * (1 - tcdf(abs(t_stat), ntime - 2));
    
    temp_T = nan(nlon*nlat, 1); temp_T(valid_idx) = trend_valid; Map_Trend(:,:,n) = reshape(temp_T, nlon, nlat)';
    temp_S = nan(nlon*nlat, 1); temp_S(valid_idx) = (p_val < 0.05); Map_Sig(:,:,n) = reshape(temp_S, nlon, nlat)';
    temp_A = nan(nlon*nlat, 1); temp_A(valid_idx) = amp_valid; Map_Amp(:,:,n) = reshape(temp_A, nlon, nlat)';
    
    Box_Trend_Ens{n} = trend_valid(:);
    Box_Amp_Ens{n}   = amp_valid(:);
    
    % --- 2. 处理独立数据集 (生成红蓝绿彩色溯源散点) ---
    for d = 1:3
        slice_raw = squeeze(TSLA_All(:,:,:,n,d));
        Y_raw = reshape(slice_raw, nlon*nlat, ntime)';
        Y_raw_valid = Y_raw(:, valid_idx);
        
        B_raw = X_mat \ Y_raw_valid;
        Box_Trend_Raw{d, n} = B_raw(1,:)' * 1000;
        Resid_raw = Y_raw_valid - X_mat * B_raw;
        Box_Amp_Raw{d, n} = std(Resid_raw, 0, 1)' * 1000;
    end
end

%% 3. 配色定义
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1 1 1; 1 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
Cmap_Trend = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

Amp_Points = [1 1 1; 1 0.9 0.5; 0.95 0.5 0.15; 0.8 0.1 0.15; 0.4 0 0]; 
Amp_Locs = [0, 0.2, 0.5, 0.8, 1]; 
Cmap_Amp = [interp1(Amp_Locs, Amp_Points(:,1), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,2), xx, 'pchip')', interp1(Amp_Locs, Amp_Points(:,3), xx, 'pchip')'];

BoxBgColor = [0.92 0.92 0.92]; % 集合大箱子用极其优雅的亮灰色
colors_raw = {
    [0.902, 0.294, 0.208], % EN4 (樱桃红)
    [0.302, 0.455, 0.690], % IAP (钴蓝)
    [0.000, 0.627, 0.408]  % Ishii (薄荷绿)
};

%% 4. 自动绘制 4 张图表
fprintf('>> 开始绘制顶刊图表...\n');

Draw_Spatial_Map(Map_Trend, Map_Sig, Lon_Grid, Lat_Grid, Cmap_Trend, sprintf('Ensemble Mean Trend Spatial Distribution (%s)', state), 'mm/yr', true);
Draw_Spatial_Map(Map_Amp, [], Lon_Grid, Lat_Grid, Cmap_Amp, sprintf('Ensemble Mean Amplitude Spatial Distribution (%s)', state), 'mm', false);

% ================= Figure 2 & 4: 统一X轴多轴箱线图 =================
Draw_Unified_MultiAxis_Boxplot(Box_Trend_Ens, Box_Trend_Raw, datasets, colors_raw, BoxBgColor, sprintf('Ensemble Trend Convergence & Sources (%s)', state), 'Trend (mm/yr)', true);
Draw_Unified_MultiAxis_Boxplot(Box_Amp_Ens, Box_Amp_Raw, datasets, colors_raw, BoxBgColor, sprintf('Ensemble Amplitude Convergence & Sources (%s)', state), 'Amplitude (mm)', false);

fprintf('======================================================\n');
fprintf('🎉 恭喜！统一X轴溯源全流程分析完毕，图表已全部生成！\n');
end

%% =========================================================================
%% 内部函数 1: 顶刊地图绘制模块
function Draw_Spatial_Map(Data_Map, Sig_Map, Lon_Grid, Lat_Grid, cmap, fig_title, unit_str, is_divergent)
    figure('Position', [50+rand*50, 50+rand*50, 1100, 1300], 'Color', 'w', 'Name', fig_title);
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:8
        nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-80 80]);
        data = Data_Map(:,:,i); m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
        
        if ~isempty(Sig_Map)
            sig = Sig_Map(:,:,i); [row, ~] = find(sig(1:3:end, 1:3:end) == 1);
            if ~isempty(row)
                lon_s = Lon_Grid(1:3:end, 1:3:end); lat_s = Lat_Grid(1:3:end, 1:3:end);
                m_plot(lon_s(sig(1:3:end, 1:3:end)==1), lat_s(sig(1:3:end, 1:3:end)==1), 'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
            end
        end
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8, 'XaxisLocation', 'bottom'); 
        colormap(gca, cmap);
        
        raw_max = prctile(abs(data(:)), 99.5); if isnan(raw_max) || raw_max==0, raw_max=1e-5; end
        exp_n = floor(log10(raw_max)); frac = raw_max / 10^exp_n;
        if frac <= 2, b = 2; elseif frac <= 5, b = 5; else, b = 10; end
        lim = b * 10^exp_n;
        
        if is_divergent, caxis([-lim, lim]); else, caxis([0, lim]); end
        title(sprintf('Term %d', i), 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location','eastoutside'); cb.FontSize = 8; title(cb, unit_str, 'FontSize', 8, 'FontWeight', 'bold');
    end
    sgtitle(fig_title, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
end

%% 内部绘图核心函数 2: 顶刊双 Y 轴箱线图模块 - 与 Plot_Combined 脚本风格一致
function Draw_Unified_MultiAxis_Boxplot(Ens_Cell, Raw_Cell, datasets, colors_raw, BoxBgColor, fig_title, y_label, show_zero_band)
    figure('Position', [50, 50, 1400, 550], 'Color', 'w', 'Name', fig_title);
    
    % 顶刊配色
    BoxWidth = 0.25;     % 减小箱体宽度，为三个数据集留出空间
    ScatterWidth = 0.2;  
    ScatterAlpha = 0.4;  
    ScatterSize  = 10;
    
    % 标准态：显示全部八阶，三个y轴
    % 调整位置，为每个数据集留出空间
    Pos_G1 = [1.2, 2.2, 3.2]; Pos_G2 = [4.4, 5.4]; Pos_G3 = [6.6, 7.6, 8.6];
    X_Limit = [0.4, 9.2];
    MainPos = [0.06, 0.15, 0.86, 0.75]; 
    
    % --- 智能极值提取 (结合集合平均和散点原始极值，保证不过曝) ---
    all_d1 = []; all_d2 = []; all_d3 = [];
    for k=1:3, all_d1 = [all_d1; Ens_Cell{k}]; for d=1:3, all_d1 = [all_d1; Raw_Cell{d,k}]; end; end
    for k=4:5, all_d2 = [all_d2; Ens_Cell{k}]; for d=1:3, all_d2 = [all_d2; Raw_Cell{d,k}]; end; end
    for k=6:8, all_d3 = [all_d3; Ens_Cell{k}]; for d=1:3, all_d3 = [all_d3; Raw_Cell{d,k}]; end; end
    
    if ~show_zero_band % 振幅必须取绝对值求边界
        all_d1 = abs(all_d1); all_d2 = abs(all_d2); all_d3 = abs(all_d3);
    end
    
    max_v1 = max(abs(prctile(all_d1, [0.1 99.9]))); if isnan(max_v1)||max_v1==0, max_v1=1; end; Limit1 = max_v1 * 1.1;
    max_v2 = max(abs(prctile(all_d2, [0.1 99.9]))); if isnan(max_v2)||max_v2==0, max_v2=1; end; Limit2 = max_v2 * 1.2;
    max_v3 = max(abs(prctile(all_d3, [0.1 99.9]))); if isnan(max_v3)||max_v3==0, max_v3=1; end; Limit3 = max_v3 * 1.2;

    % --- Layer 1 (Red) ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for k = 1:3
        data = Ens_Cell{k}; if isempty(data), continue; end
        % 振幅数据确保为正
        if ~show_zero_band
            data = abs(data);
        end
        
        % 为每个数据集绘制箱体和散点
        for i = 1:3
            dataset_data = Raw_Cell{i, k}; if isempty(dataset_data), continue; end
            % 振幅数据确保为正
            if ~show_zero_band
                dataset_data = abs(dataset_data);
            end
            % 计算每个数据集的x位置，在同一根x轴上并排显示
            x_pos = Pos_G1(k) - 0.3 + (i-1)*0.3;
            % 绘制数据集箱体和散点
            draw_jitter_unit(ax1, x_pos, dataset_data, BoxWidth, ScatterWidth, ScatterSize, colors_raw{i}, BoxBgColor, ScatterAlpha);
        end
    end
    % 智能判断范围：趋势用 [-Max, Max]，振幅用 [0, Max]
    if show_zero_band
        ylim(ax1, [-Limit1, Limit1]); 
        ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
    else
        ylim(ax1, [0, Limit1]); 
        ylabel(ax1, 'Dominant Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    end
    set(ax1, 'XLim', X_Limit, 'YColor', 'k', 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- Layer 2 (Blue) ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        data = Ens_Cell{k+3}; if isempty(data), continue; end
        % 振幅数据确保为正
        if ~show_zero_band
            data = abs(data);
        end
        
        % 为每个数据集绘制箱体和散点
        for i = 1:3
            dataset_data = Raw_Cell{i, k+3}; if isempty(dataset_data), continue; end
            % 振幅数据确保为正
            if ~show_zero_band
                dataset_data = abs(dataset_data);
            end
            % 计算每个数据集的x位置，在同一根x轴上并排显示
            x_pos = Pos_G2(k) - 0.3 + (i-1)*0.3;
            % 绘制数据集箱体和散点
            draw_jitter_unit(ax2, x_pos, dataset_data, BoxWidth, ScatterWidth, ScatterSize, colors_raw{i}, BoxBgColor, ScatterAlpha);
        end
    end
    if show_zero_band
        ylim(ax2, [-Limit2, Limit2]); 
        ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
    else
        ylim(ax2, [0, Limit2]); 
        ylabel(ax2, 'Secondary Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    end
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', 'k', 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- Layer 3 (Green) ---
    Offset_Ratio = 0.08; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = X_Limit * Scale_Factor; 
    
    for k = 1:3
        data = Ens_Cell{k+5}; if isempty(data), continue; end
        % 振幅数据确保为正
        if ~show_zero_band
            data = abs(data);
        end
        
        % 为每个数据集绘制箱体和散点
        for i = 1:3
            dataset_data = Raw_Cell{i, k+5}; if isempty(dataset_data), continue; end
            % 振幅数据确保为正
            if ~show_zero_band
                dataset_data = abs(dataset_data);
            end
            % 计算每个数据集的x位置，在同一根x轴上并排显示
            x_pos = Pos_G3(k) - 0.3 + (i-1)*0.3;
            % 绘制数据集箱体和散点
            draw_jitter_unit(ax3, x_pos, dataset_data, BoxWidth, ScatterWidth, ScatterSize, colors_raw{i}, BoxBgColor, ScatterAlpha);
        end
    end
    if show_zero_band
        ylim(ax3, [-Limit3, Limit3]);
        ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 16);
    else
        ylim(ax3, [0, Limit3]);
        ylabel(ax3, 'High-Order Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    end
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', 'k', 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % --- 装饰与分割线 ---
    yline(ax1, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.2); 
    xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    
    % 设置 X 轴刻度和标签，显示每个阶数的位置
    All_Ticks = [1.2, 2.2, 3.2, 4.4, 5.4, 6.6, 7.6, 8.6];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    title(ax1, fig_title, 'FontSize', 18, 'FontWeight', 'bold');
    
    % === 图例 (右上角纯白无边框紧凑版) ===
    h_en4  = scatter(ax1, NaN, NaN, ScatterSize, colors_raw{1}, 'filled', 'MarkerFaceAlpha', 0.6);
    h_iap  = scatter(ax1, NaN, NaN, ScatterSize, colors_raw{2}, 'filled', 'MarkerFaceAlpha', 0.6);
    h_ishii = scatter(ax1, NaN, NaN, ScatterSize, colors_raw{3}, 'filled', 'MarkerFaceAlpha', 0.6);
    
    legend(ax1, [h_en4, h_iap, h_ishii], ...
           {'EN4', 'IAP', 'Ishii'}, ...
           'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');
    
    % 添加数据集说明文本
    text(ax1, 0.02, 0.95, '每个阶数显示三个数据集:', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
    text(ax1, 0.02, 0.90, '● EN4 (红色) | ● IAP (蓝色) | ● Ishii (绿色)', 'Units', 'normalized', 'FontSize', 11);

    % === 交互优化 ===
    Ratio2 = Limit2 / Limit1; Ratio3 = Limit3 / Limit1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor, show_zero_band);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor, show_zero_band);
end

%% 内部函数: 绘制抖动箱线图单元
function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    
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
    
    % 绘制箱体 - 使用统一的灰色箱体
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], [0.7 0.7 0.7], ...
        'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    % 绘制中位线和均值点
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

%% 内部辅助: 同步多轴缩放
function SyncAllAxes(ax1, ax2, ax3, R2, R3, Factor, show_zero_band)
    try
        if show_zero_band
            L = max(abs(ax1.YLim));
            ax1.YLim = [-L, L]; ax2.YLim = [-L*R2, L*R2]; ax3.YLim = [-L*R3, L*R3];
        else
            L = max(ax1.YLim);
            ax1.YLim = [0, L]; ax2.YLim = [0, L*R2]; ax3.YLim = [0, L*R3];
        end
        ax3.XLim = ax1.XLim * Factor;
    catch
    end
end