function Plot_Amplitude_Combined(dataset_name, state, target_comp)
% 功能: 整合绘制振幅空间分布图和箱线图 (顶刊标准)
% 特性: 动态自适应阶数，完美兼容 8阶(TSLA) 与 单阶(SSLA) 数据
% 输入参数:
%   dataset_name: 数据集名称 ('EN4', 'IAP', 'Ishii')
%   state: 状态 ('StdRef', 'Average')
%   target_comp: 要绘制的变量名称 ('TSLA', 'HSLA', 'SSLA', 或 'Auto') [可选，默认 'Auto']

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

if nargin < 3
    target_comp = 'Auto';
end

%% 2. 准备数据 (根据输入选择对应变量，或兼容默认的优先逻辑)
if (strcmpi(target_comp, 'TSLA') || strcmpi(target_comp, 'Auto')) && exist('amp_TSLA', 'var') && ~isempty(amp_TSLA)
    amplitude_data = amp_TSLA;
    comp_name = 'TSLA (Thermosteric)';
    fprintf('>> ✅ 选择绘制 TSLA 的振幅图。\n');
elseif (strcmpi(target_comp, 'HSLA') || strcmpi(target_comp, 'Auto')) && exist('amp_HSLA', 'var') && ~isempty(amp_HSLA)
    amplitude_data = amp_HSLA;
    comp_name = 'HSLA (Halosteric)';
    fprintf('>> ✅ 选择绘制 HSLA 的振幅图。\n');
elseif (strcmpi(target_comp, 'SSLA') || strcmpi(target_comp, 'Auto')) && exist('amp_SSLA', 'var') && ~isempty(amp_SSLA)
    amplitude_data = amp_SSLA;
    comp_name = 'SSLA (Total Steric)';
    fprintf('>> ✅ 选择绘制 SSLA 的振幅图。\n');
else
    error('❌ 缺少请求的振幅变量 (%s) 或无法找到任何有效变量！', target_comp);
end

% ！！！智能单位换算 (m -> mm) ！！！
if max(abs(amplitude_data(:))) < 0.2
    fprintf('>> 🔄 检测到单位可能为 [米(m)]，自动乘以 1000 转换为 [毫米(mm)]...\n');
    amplitude_data = amplitude_data * 1000;
end

% 经纬度处理
if ~exist('Lon', 'var') && exist('target_lon', 'var'), Lon = target_lon; end
if ~exist('Lat', 'var') && exist('target_lat', 'var'), Lat = target_lat; end
if ~exist('Lon', 'var'), error('❌ 未找到经度数据'); end
if ~exist('Lat', 'var'), error('❌ 未找到纬度数据'); end

% 确保Lon和Lat是向量
if ~isvector(Lon), Lon = Lon(1,:); end
if ~isvector(Lat), Lat = Lat(:,1); end

% ⚠️ 第一步：先转置使得矩阵符合地图绘制要求 [Lat x Lon]
if size(amplitude_data, 1) ~= length(Lat)
    temp_a = nan(length(Lat), length(Lon), size(amplitude_data, 3));
    for i = 1:size(amplitude_data, 3)
        temp_a(:,:,i) = amplitude_data(:,:,i)';
    end
    amplitude_data = temp_a;
end

% ⚠️ 第二步：再进行经度 [0, 360] -> [-180, 180] 转换与切片拼接
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    for i = 1:size(amplitude_data, 3)
        amplitude_data(:,:,i) = amplitude_data(:, sort_idx, i);
    end
end

[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

%% 3. 绘制空间分布图 (顶刊 振幅热力 色棒)
fprintf('>> 绘制振幅空间分布图...\n');

max_terms = size(amplitude_data, 3);
if max_terms > 8, max_terms = 8; end % 最多画8阶

figure('Position', [100, 50, 1200, 1400], 'Color', 'w', 'Name', 'Spatial Maps');
t = tiledlayout(ceil(max_terms/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

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

sgtitle(sprintf('%s (%s) 1-%d Terms Amplitude Distribution [%s]', dataset_name, state, max_terms, comp_name), 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% 4. 准备箱线图数据
fprintf('>> 准备箱线图数据...\n');
mask = ~isnan(amplitude_data(:,:,1));
[rows, cols] = find(mask);
num_valid = length(rows);
Data_Amplitude = nan(num_valid, max_terms);

for n = 1:max_terms
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

%% 5. 绘制箱线图 (🚨 自适应阶数保护机制)
fprintf('>> 绘制振幅箱线图...\n');
figure('Position', [50, 50, 1200, 550], 'Color', 'w', 'Name', 'Amplitude Boxplot');

% 顶刊配色
C1 = [0.902, 0.294, 0.208]; % 樱桃红
C2 = [0.302, 0.455, 0.690]; % 钴蓝
C3 = [0.000, 0.627, 0.408]; % 薄荷绿
BoxBgColor = [0.96 0.96 0.96]; % 极浅高级灰底色

BoxWidth = 0.55;     
ScatterWidth = 0.4;  
ScatterAlpha = 0.4;  
ScatterSize  = 12;

X_Limit = [0.4, 9.0];
MainPos = [0.07, 0.10, 0.80, 0.73]; % 进一步降低绘图区高度，给上方留足空间

if max_terms == 8
    % ==========================================
    % 当存在 8 阶数据时，执行完美堆叠箱线图逻辑
    % ==========================================
    if strcmp(state, 'Average')
        Pos_G1 = [1]; Pos_G2 = [2, 3]; Pos_G3 = [4.2, 5.2]; Pos_G4 = [6.4, 7.4, 8.4]; 
        
        ax1 = axes('Position', MainPos); hold(ax1, 'on');
        draw_jitter_unit(ax1, Pos_G1(1), Data_Amplitude(:,1), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
        max_v1 = max(Data_Amplitude(:,1)); 
        if isnan(max_v1)||max_v1==0, max_v1=1; end
        ylim(ax1, [0, max_v1 * 1.1]); 
        set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        ylabel(ax1, 'Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
        
        ax2 = axes('Position', MainPos); hold(ax2, 'on');
        for k = 1:2
            draw_jitter_unit(ax2, Pos_G2(k), Data_Amplitude(:,k+1), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
        end
        max_v2 = max(Data_Amplitude(:,2:3), [], 'all'); 
        if isnan(max_v2)||max_v2==0, max_v2=1; end
        ylim(ax2, [0, max_v2 * 1.2]); 
        set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        
        Offset_Ratio = 0.06; 
        Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
        ax3 = axes('Position', Pos3); hold(ax3, 'on');
        Scale_Factor = Pos3(3) / MainPos(3); 
        Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factor]; 
        
        for k = 1:2
            draw_jitter_unit(ax3, Pos_G3(k), Data_Amplitude(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
        end
        max_v3 = max(Data_Amplitude(:,4:5), [], 'all'); 
        if isnan(max_v3)||max_v3==0, max_v3=1; end
        ylim(ax3, [0, max_v3 * 1.2]);
        set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        
        Offset_Ratio2 = 0.12; 
        Pos4 = MainPos; Pos4(3) = MainPos(3) * (1 + Offset_Ratio2);
        ax4 = axes('Position', Pos4); hold(ax4, 'on');
        Scale_Factor2 = Pos4(3) / MainPos(3); 
        Real_X_Limit2 = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factor2]; 
        
        for k = 1:3
            draw_jitter_unit(ax4, Pos_G4(k), Data_Amplitude(:,k+5), BoxWidth, ScatterWidth, ScatterSize, [0.5 0 0.5], BoxBgColor, ScatterAlpha);
        end
        max_v4 = max(Data_Amplitude(:,6:8), [], 'all'); 
        if isnan(max_v4)||max_v4==0, max_v4=1; end
        ylim(ax4, [0, max_v4 * 1.2]);
        set(ax4, 'Color', 'none', 'XLim', Real_X_Limit2, 'YAxisLocation', 'right', 'YColor', [0.5 0 0.5], 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        
        xline(ax1, 1.5, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
        xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
        xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
        yline(ax1, 0, '-', 'Color', [0 0 0], 'LineWidth', 1.5); % 绘制强力底部 0 轴
        
        All_Ticks = [Pos_G1, Pos_G2, Pos_G3, Pos_G4]; 
        c1s = sprintf('\\color[rgb]{%f,%f,%f}', C1);
        c2s = sprintf('\\color[rgb]{%f,%f,%f}', C2);
        c3s = sprintf('\\color[rgb]{%f,%f,%f}', C3);
        c4s = sprintf('\\color[rgb]{%f,%f,%f}', 0.5, 0, 0.5);
        All_Labels = {[c1s 'T1'], [c2s 'T2'], [c2s 'T3'], [c3s 'T4'], [c3s 'T5'], [c4s 'T6'], [c4s 'T7'], [c4s 'T8']};
        set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold', 'TickLabelInterpreter', 'tex');
        
        Ratio2 = max_v2 / max_v1; Ratio3 = max_v3 / max_v1; Ratio4 = max_v4 / max_v1;
        linkaxes([ax1, ax2], 'x'); 
        hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
        hZoom.ActionPostCallback = @(s,e) SyncAllAxes4(ax1, ax2, ax3, ax4, Ratio2, Ratio3, Ratio4, Scale_Factor, Scale_Factor2);
        hPan = pan(gcf); hPan.Motion = 'both';
        hPan.ActionPostCallback = @(s,e) SyncAllAxes4(ax1, ax2, ax3, ax4, Ratio2, Ratio3, Ratio4, Scale_Factor, Scale_Factor2);
    else
        % 标准态划分
        Pos_G1 = [1, 2, 3]; Pos_G2 = [4.2, 5.2]; Pos_G3 = [6.4, 7.4, 8.4]; 
        ax1 = axes('Position', MainPos); hold(ax1, 'on');
        for k = 1:3
            draw_jitter_unit(ax1, Pos_G1(k), Data_Amplitude(:,k), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
        end
        max_v1 = max(Data_Amplitude(:,1:3), [], 'all'); 
        if isnan(max_v1)||max_v1==0, max_v1=1; end
        ylim(ax1, [0, max_v1 * 1.1]); 
        set(ax1, 'XLim', X_Limit, 'YColor', C1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        ylabel(ax1, 'Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
        
        ax2 = axes('Position', MainPos); hold(ax2, 'on');
        for k = 1:2
            draw_jitter_unit(ax2, Pos_G2(k), Data_Amplitude(:,k+3), BoxWidth, ScatterWidth, ScatterSize, C2, BoxBgColor, ScatterAlpha);
        end
        max_v2 = max(Data_Amplitude(:,4:5), [], 'all'); 
        if isnan(max_v2)||max_v2==0, max_v2=1; end
        ylim(ax2, [0, max_v2 * 1.2]); 
        set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        
        Offset_Ratio = 0.04; 
        Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
        ax3 = axes('Position', Pos3); hold(ax3, 'on');
        Scale_Factor = Pos3(3) / MainPos(3); 
        Real_X_Limit = [X_Limit(1), X_Limit(1) + (X_Limit(2) - X_Limit(1)) * Scale_Factor]; 
        
        for k = 1:3
            draw_jitter_unit(ax3, Pos_G3(k), Data_Amplitude(:,k+5), BoxWidth, ScatterWidth, ScatterSize, C3, BoxBgColor, ScatterAlpha);
        end
        max_v3 = max(Data_Amplitude(:,6:8), [], 'all'); 
        if isnan(max_v3)||max_v3==0, max_v3=1; end
        ylim(ax3, [0, max_v3 * 1.2]);
        set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        
        xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5); 
        xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
        yline(ax1, 0, '-', 'Color', [0 0 0], 'LineWidth', 1.5); % 绘制强力底部 0 轴
        
        All_Ticks = [Pos_G1, Pos_G2, Pos_G3]; 
        c1s = sprintf('\\color[rgb]{%f,%f,%f}', C1);
        c2s = sprintf('\\color[rgb]{%f,%f,%f}', C2);
        c3s = sprintf('\\color[rgb]{%f,%f,%f}', C3);
        All_Labels = {[c1s 'T1'], [c1s 'T2'], [c1s 'T3'], [c2s 'T4'], [c2s 'T5'], [c3s 'T6'], [c3s 'T7'], [c3s 'T8']};
        set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold', 'TickLabelInterpreter', 'tex');
        
        Ratio2 = max_v2 / max_v1; Ratio3 = max_v3 / max_v1;
        linkaxes([ax1, ax2], 'x'); 
        hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
        hZoom.ActionPostCallback = @(s,e) SyncThreeAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
        hPan = pan(gcf); hPan.Motion = 'both';
        hPan.ActionPostCallback = @(s,e) SyncThreeAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    end

elseif max_terms == 1
    % ==========================================
    % 防撞墙兜底：当只加载了 SSLA (1阶) 时的单箱图画法
    % ==========================================
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    draw_jitter_unit(ax1, 1, Data_Amplitude(:,1), BoxWidth, ScatterWidth, ScatterSize, [0.7 0.1 0.1], BoxBgColor, ScatterAlpha);
    max_v1 = max(Data_Amplitude(:,1)); 
    if isnan(max_v1)||max_v1==0, max_v1=1; end
    ylim(ax1, [0, max_v1 * 1.1]); 
    set(ax1, 'XLim', [0.4, 1.6], 'YColor', C1, 'Box', 'on', 'XTick', 1, 'XTickLabel', {sprintf('\\color[rgb]{%f,%f,%f}Total Amplitude', C1)}, 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'TickLabelInterpreter', 'tex');
    ylabel(ax1, 'Amplitude (mm)', 'FontWeight', 'bold', 'FontSize', 16);
end

% ==========================================================
% 统一通用图例与标题 (使用隐形大坐标轴完美挂载标题与图例，杜绝遮挡)
% ==========================================================
xlabel(ax1, sprintf('%s Components', strtok(comp_name, ' ')), 'FontSize', 16, 'FontWeight', 'bold');

% 创建一个充满全图顶部的透明 axes，专门用来放 Title
ax_title = axes('Position', [0, 0.90, 1, 0.1], 'Visible', 'off');
text(ax_title, 0.5, 0.8, sprintf('%s (%s) %s Amplitude Distribution', dataset_name, state, comp_name), ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');

h_med     = plot(ax1, NaN, NaN, 'r-', 'LineWidth', 1.5);
h_mean    = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h_whisker = plot(ax1, NaN, NaN, 'k-', 'LineWidth', 1.5);
h_out     = scatter(ax1, NaN, NaN, ScatterSize, C1, 'filled', 'MarkerFaceAlpha', 0.8);

outlier_str = sprintf('Outliers-%s', dataset_name);

% ==========================================================
% 自定义图例 (完美对齐 X 轴 && 带红框，散落分布)
% ==========================================================
% 画出与主图同宽的自定义坐标系充当图例框 (去除边框)
ax_lgd = axes('Position', [MainPos(1), 0.88, MainPos(3), 0.025], ...
    'XLim', [0, 100], 'YLim', [0, 1], ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none', ...
    'XColor', 'none', 'YColor', 'none');
hold(ax_lgd, 'on');

% 第 1 项: Median line (红线)
plot(ax_lgd, [6, 11], [0.5, 0.5], 'r-', 'LineWidth', 2);
text(ax_lgd, 12, 0.5, 'Median line', 'FontSize', 14, 'FontName', 'Times New Roman', 'VerticalAlignment', 'middle');

% 第 2 项: Mean value (黑点)
plot(ax_lgd, 31, 0.5, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
text(ax_lgd, 33, 0.5, 'Mean value', 'FontSize', 14, 'FontName', 'Times New Roman', 'VerticalAlignment', 'middle');

% 第 3 项: 2.5% ~ 97.5% (黑线)
plot(ax_lgd, [53, 58], [0.5, 0.5], 'k-', 'LineWidth', 1.5);
text(ax_lgd, 59, 0.5, '2.5% \sim 97.5%', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', 'Times New Roman', 'VerticalAlignment', 'middle');

% 第 4 项: Outliers (彩色点，根据当前数据集)
scatter(ax_lgd, 80, 0.5, 40, C1, 'filled', 'MarkerFaceAlpha', 0.8);
text(ax_lgd, 82, 0.5, outlier_str, 'FontSize', 14, 'FontName', 'Times New Roman', 'VerticalAlignment', 'middle');

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
    
    % 绘制中位线、均值点
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
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

function SyncAllAxes4(ax1, ax2, ax3, ax4, R2, R3, R4, Factor1, Factor2)
    try
        L = max(ax1.YLim);
        ax1.YLim = [0, L];
        ax2.YLim = [0, L*R2];
        ax3.YLim = [0, L*R3];
        ax4.YLim = [0, L*R4];
        ax3.XLim = [ax1.XLim(1), ax1.XLim(1) + (ax1.XLim(2) - ax1.XLim(1)) * Factor1];
        ax4.XLim = [ax1.XLim(1), ax1.XLim(1) + (ax1.XLim(2) - ax1.XLim(1)) * Factor2];
    catch ME
        disp(['SyncAllAxes4 错误: ', ME.message]);
    end
end

function SyncThreeAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        L = max(ax1.YLim);
        ax1.YLim = [0, L];
        ax2.YLim = [0, L*R2];
        ax3.YLim = [0, L*R3];
        ax3.XLim = [ax1.XLim(1), ax1.XLim(1) + (ax1.XLim(2) - ax1.XLim(1)) * Factor];
    catch ME
        disp(['SyncThreeAxes 错误: ', ME.message]);
    end
end