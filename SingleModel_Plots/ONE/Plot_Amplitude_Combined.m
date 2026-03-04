function Plot_Amplitude_Combined(dataset_name, state)
% 功能: 整合绘制振幅空间分布图和箱线图
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')

clc; close all;
fprintf('========== 整合绘制振幅空间分布图和箱线图 ==========\n');

%% 1. 加载振幅结果文件
amplitude_dir = 'D:\work\Task_Convergence\Amplitude_Results';
amplitude_file = fullfile(amplitude_dir, sprintf('%s_%s_Amplitude.mat', dataset_name, state));

if ~exist(amplitude_file, 'file')
    error('❌ 找不到振幅结果文件: %s', amplitude_file);
end

fprintf('>> 加载振幅数据: %s ...\n', amplitude_file);
load(amplitude_file);

%% 2. 准备数据
% 提取振幅数据 (优先使用 HSLA 以显示真实信号，若无则使用 TSLA)
if exist('amp_HSLA', 'var') && ~isempty(amp_HSLA)
    amplitude_data = amp_HSLA;
    fprintf('>> 检测到 HSLA (总比容)，将绘制 HSLA 振幅。\n');
else
    amplitude_data = amp_TSLA;
    fprintf('>> 未检测到 HSLA，将绘制纯温度 TSLA 振幅。\n');
end

% 经纬度处理
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end

if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% 坐标转换 [0, 360] -> [-180, 180]
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    for i = 1:size(amplitude_data, 3)
        amplitude_data(:,:,i) = amplitude_data(:, sort_idx, i);
    end
end

[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

% 转置使得画图时不报错
if size(amplitude_data, 1) ~= length(Lat)
    for i = 1:size(amplitude_data, 3)
        temp_a(:,:,i) = amplitude_data(:,:,i)';
    end
    amplitude_data = temp_a;
end

%% 3. 绘制空间分布图 (顶刊 振幅热力 色棒)
fprintf('>> 绘制振幅空间分布图...\n');

% 确定显示的阶数
if strcmp(state, 'Average')
    max_terms = 8;
    figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Spatial Maps');
    t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    fprintf('>> 检测到平均态，绘制全部八阶振幅...\n');
else
    max_terms = 8;
    figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Spatial Maps');
    t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    fprintf('>> 检测到标准态，绘制全部八阶振幅...\n');
end

% 顶刊振幅专属色盘：白 -> 淡黄 -> 亮橘 -> 深红 
Amp_Points = [ 
    1.00, 1.00, 1.00;  % 0.0: 纯白 (背景/极弱区) 
    1.00, 0.90, 0.50;  % 0.2: 暖黄 (弱振幅) 
    0.95, 0.50, 0.15;  % 0.5: 亮橙 (中等振幅) 
    0.80, 0.10, 0.15;  % 0.8: 强红 (强振幅) 
    0.40, 0.00, 0.00   % 1.0: 暗红 (极值核心区) 
]; 
Amp_Locs = [0, 0.2, 0.5, 0.8, 1]; 
xx = linspace(0, 1, 256); 
MyAmpCmap = [interp1(Amp_Locs, Amp_Points(:,1), xx, 'pchip')', ... 
             interp1(Amp_Locs, Amp_Points(:,2), xx, 'pchip')', ... 
             interp1(Amp_Locs, Amp_Points(:,3), xx, 'pchip')'];

for term_idx = 1:max_terms
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    data = abs(amplitude_data(:,:,term_idx)); % 确保振幅为正
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
    
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    
    colormap(gca, MyAmpCmap);
    
    raw_max = prctile(data(:), 99.5); % 振幅用 99.5% 截断更佳
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    
    % ========================================================
    % 🚀 核心修复 1：振幅是正值，必须使用 [0, Max]！
    % ========================================================
    caxis([0, final_limit]); 
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    cb = colorbar('Location', 'eastoutside'); cb.FontSize = 8;
    title(cb, 'mm', 'FontSize', 8, 'FontWeight', 'bold'); 
    cb.Ticks = linspace(0, final_limit, 5); % Ticks 同步修改为正数
end

sgtitle(sprintf('%s (%s) 1-%d Terms Amplitude Distribution', dataset_name, state, max_terms), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 4. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(amplitude_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Amplitude = nan(num_valid, 8);

for n = 1:8
    data_slice = amplitude_data(:,:,n);
    vec_a = nan(num_valid, 1);
    for k = 1:num_valid
        vec_a(k) = abs(data_slice(rows(k), cols(k)));
    end
    Data_Amplitude(:,n) = vec_a;
end

if max(Data_Amplitude(:)) > 500
     Data_Amplitude = Data_Amplitude / 1000;
end

%% 5. 绘制箱线图
fprintf('>> 绘制振幅箱线图...\n');
figure('Position', [50, 50, 1100, 500], 'Color', 'w', 'Name', 'Amplitude Boxplot');

% 顶刊配色
C1 = [0.902, 0.294, 0.208]; % 樱桃红
C2 = [0.302, 0.455, 0.690]; % 钴蓝
C3 = [0.000, 0.627, 0.408]; % 薄荷绿
BoxBgColor = [0.96 0.96 0.96]; % 极浅高级灰底色

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

if strcmp(state, 'Average')
    % 平均态：显示全部八阶，三个y轴
    Pos_G1 = [1, 2, 3]; Pos_G2 = [4.5, 5.5]; Pos_G3 = [6.8, 7.8, 8.8];
    X_Limit = [0.4, 9.5];
    MainPos = [0.08, 0.15, 0.82, 0.75]; 
    
    % --- Layer 1 (Red) ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for k = 1:3
        draw_jitter_unit(ax1, Pos_G1(k), Data_Amplitude(:,k), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
    end
    max_v1 = max(Data_Amplitude(:,1:3), [], 'all'); 
    if isnan(max_v1)||max_v1==0, max_v1=1; end
    Limit1 = max_v1 * 1.1; 
    ylim(ax1, [0, Limit1]); 
    set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax1, 'Dominant Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 2 (Blue) ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        draw_jitter_unit(ax2, Pos_G2(k), Data_Amplitude(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    max_v2 = max(Data_Amplitude(:,4:5), [], 'all'); 
    if isnan(max_v2)||max_v2==0, max_v2=1; end
    Limit2 = max_v2 * 1.2;
    ylim(ax2, [0, Limit2]); 
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax2, 'Secondary Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 3 (Green) ---
    Offset_Ratio = 0.08; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factor]; 
    
    for k = 1:3
        draw_jitter_unit(ax3, Pos_G3(k), Data_Amplitude(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
    end
    max_v3 = max(Data_Amplitude(:,6:8), [], 'all'); 
    if isnan(max_v3)||max_v3==0, max_v3=1; end
    Limit3 = max_v3 * 1.2;
    ylim(ax3, [0, Limit3]);
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax3, 'High-Order Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- 装饰与分割线 ---
    xline(ax1, 3.75, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 6.1, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    yline(ax1, 0, '-', 'Color', [0 0 0], 'LineWidth', 1.5); % 绘制强力底部 0 轴
    
    All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    title(ax1, sprintf('%s (%s) Terms Amplitude Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');
    
    % === 交互优化 ===
    Ratio2 = Limit2 / Limit1; Ratio3 = Limit3 / Limit1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    
else
    % 标准态：显示全部八阶，三个y轴
    Pos_G1 = [1, 2, 3]; Pos_G2 = [4.5, 5.5]; Pos_G3 = [6.8, 7.8, 8.8];
    X_Limit = [0.4, 9.5];
    MainPos = [0.08, 0.15, 0.82, 0.75]; 
    
    % --- Layer 1 (Red) ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for k = 1:3
        draw_jitter_unit(ax1, Pos_G1(k), Data_Amplitude(:,k), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
    end
    max_v1 = max(Data_Amplitude(:,1:3), [], 'all'); 
    if isnan(max_v1)||max_v1==0, max_v1=1; end
    Limit1 = max_v1 * 1.1; 
    ylim(ax1, [0, Limit1]); 
    set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax1, 'Dominant Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 2 (Blue) ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        draw_jitter_unit(ax2, Pos_G2(k), Data_Amplitude(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
    end
    max_v2 = max(Data_Amplitude(:,4:5), [], 'all'); 
    if isnan(max_v2)||max_v2==0, max_v2=1; end
    Limit2 = max_v2 * 1.2;
    ylim(ax2, [0, Limit2]); 
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax2, 'Secondary Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- Layer 3 (Green) ---
    Offset_Ratio = 0.08; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factor]; 
    
    for k = 1:3
        draw_jitter_unit(ax3, Pos_G3(k), Data_Amplitude(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
    end
    max_v3 = max(Data_Amplitude(:,6:8), [], 'all'); 
    if isnan(max_v3)||max_v3==0, max_v3=1; end
    Limit3 = max_v3 * 1.2;
    ylim(ax3, [0, Limit3]);
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax3, 'High-Order Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
    
    % --- 装饰与分割线 ---
    xline(ax1, 3.75, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 6.1, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    yline(ax1, 0, '-', 'Color', [0 0 0], 'LineWidth', 1.5); % 绘制强力底部 0 轴
    
    All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Thermosteric Components (Terms)', 'FontSize', 16, 'FontWeight', 'bold');
    title(ax1, sprintf('%s (%s) Terms Amplitude Distribution', dataset_name, state), 'FontSize', 18, 'FontWeight', 'bold');
    
    % === 交互优化 ===
    Ratio2 = Limit2 / Limit1; Ratio3 = Limit3 / Limit1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
end

% === 图例 (右上角纯白无边框紧凑版) ===
h_box  = patch(ax1, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
h_med  = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 2);
h_mean = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_max  = plot(ax1, NaN, NaN, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % 最大值点图例
h_out  = scatter(ax1, NaN, NaN, ScatterSize, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.6); % 统一灰色代表散点

legend(ax1, [h_box, h_med, h_mean, h_max, h_out], ...
       {'IQR Box (25%-75%)', 'Median', 'Mean', 'Max', 'Outliers'}, ...
       'Location', 'northeast', 'FontSize', 11, 'Box', 'on', 'Color', [1 1 1], 'EdgeColor', 'none');

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

function SyncAllAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        L = max(ax1.YLim);
        ax1.YLim = [0, L];
        ax2.YLim = [0, L*R2];
        ax3.YLim = [0, L*R3];
        ax3.XLim = [ax1.XLim(1), ax1.XLim(1) + (ax1.XLim(2) - ax1.XLim(1)) * Factor];
    catch
    end
end

function SyncTwoAxes(ax1, ax2, R2)
    try
        L = max(ax1.YLim);
        ax1.YLim = [0, L];
        ax2.YLim = [0, L*R2];
    catch
    end
end