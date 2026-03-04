function Plot_Advanced_Boxplot()
% 功能：展示四种顶刊箱线图优化方案 (支持开关控制)
% 推荐使用 Standard Reference 数据测试，效果最震撼

clear; clc; close all;

% =========================================================
% 🎛️ 优化方案控制面板 (修改这里看各自的效果)
% =========================================================
Opt_ShadowBand = true;  % 优化一：±0.01 收敛阈值阴影带
Opt_Raincloud  = true;  % 优化二：雨云图 (右侧核密度小提琴 + 左侧散点)
Opt_MeanCurve  = true;  % 优化三：均值收敛轨迹连线
Opt_SymLog     = true; % 优化四：对称对数坐标 (消除极大极小值的量级鸿沟)
% =========================================================

%% 1. 读取数据 (以 EN4 标准态为例)
dataset = 'EN4'; 
state = 'StdRef'; 
trend_file = fullfile('D:\work\Task_Convergence\Trend_Results', sprintf('%s_%s_Trends.mat', dataset, state));

if ~exist(trend_file, 'file')
    error('找不到趋势文件: %s\n请确认路径是否正确。', trend_file);
end

fprintf('>> 加载数据: %s (%s)...\n', dataset, state);
load(trend_file);

% 提取数据 (优先用 HSLA，无则用 TSLA)
if exist('trend_HSLA', 'var') && ~isempty(trend_HSLA)
    raw_data = trend_HSLA;
else
    raw_data = trend_TSLA;
end

% 展平为1D数组，剔除 NaN
Data_Terms = cell(1, 8);
for i = 1:8
    slice = raw_data(:,:,i);
    Data_Terms{i} = slice(~isnan(slice(:)));
end

%% 2. 对称对数变换函数 (伪 Symlog)
% 如果开启 Symlog，用此函数处理极端量级
C = 0.01; % 线性与对数的过渡阈值
trans = @(y) sign(y) .* log10(1 + abs(y)/C);
inv_trans = @(y_trans) sign(y_trans) .* C .* (10.^abs(y_trans) - 1);

if Opt_SymLog
    for i = 1:8
        Data_Terms{i} = trans(Data_Terms{i});
    end
end

%% 3. 开始绘图
figure('Position', [150, 150, 1000, 600], 'Color', 'w', 'Name', 'Advanced Boxplot');
ax = axes('Position', [0.1, 0.15, 0.85, 0.75]); hold on;

% 顶刊配色
Color_Main = [0.902, 0.294, 0.208]; % 樱桃红
Color_Cloud = [0.961, 0.612, 0.420]; % 浅桃色 (用于雨云密度图)

% --- 优化一：收敛阴影带 ---
if Opt_ShadowBand
    if Opt_SymLog
        y_lim_band = [trans(-0.01), trans(0.01)];
    else
        y_lim_band = [-0.01, 0.01];
    end
    % 画一条贯穿全图的浅灰色透明带
    patch([0, 9, 9, 0], [y_lim_band(1), y_lim_band(1), y_lim_band(2), y_lim_band(2)], ...
          [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    % 添加标注文字
    text(8.8, y_lim_band(2) + diff(y_lim_band)*0.5, 'Convergence Zone (\pm0.01)', ...
         'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'right', 'FontAngle', 'italic');
end

% 画 0 轴
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% --- 初始化均值记录数组 (用于优化三) ---
means_val = zeros(1, 8);

% --- 循环绘制每一阶 ---
for i = 1:8
    data = Data_Terms{i};
    means_val(i) = mean(data);
    
    % 计算统计量
    q1 = prctile(data, 25); q3 = prctile(data, 75); med = median(data);
    p2_5 = prctile(data, 2.5); p97_5 = prctile(data, 97.5);
    
    % X 轴中心位置
    center = i;
    
    if Opt_Raincloud
        % ----------------------------------------------------
        % 优化二：雨云图 (Raincloud Plot)
        % ----------------------------------------------------
        % 1. 画云 (右侧核密度估计)
        try
            [f, xi] = ksdensity(data);
            f = f / max(f) * 0.35; % 控制云朵最大宽度为 0.35
            fill(center + 0.05 + f, xi, Color_Cloud, 'FaceAlpha', 0.6, 'EdgeColor', Color_Main, 'LineWidth', 1.2);
        catch
        end
        
        % 2. 画雨 (左侧散点)
        outliers = data(data < p2_5 | data > p97_5);
        if ~isempty(outliers)
            if length(outliers) > 1000, outliers = outliers(randperm(length(outliers), 1000)); end
            x_jit = center - 0.2 + (rand(size(outliers)) - 0.5) * 0.15;
            scatter(x_jit, outliers, 10, Color_Main, 'filled', 'MarkerFaceAlpha', 0.3);
        end
        
        % 3. 画狭窄箱体 (中间)
        patch([center-0.02, center+0.02, center+0.02, center-0.02], [q1, q1, q3, q3], ...
              [0.2 0.2 0.2], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        plot([center-0.02, center+0.02], [med, med], 'w-', 'LineWidth', 1.5); % 白线为中位数
        plot([center, center], [p2_5, q1], 'k-', 'LineWidth', 1);
        plot([center, center], [q3, p97_5], 'k-', 'LineWidth', 1);
        
    else
        % ----------------------------------------------------
        % 原始样式：标准抖动箱线图
        % ----------------------------------------------------
        outliers = data(data < p2_5 | data > p97_5);
        if ~isempty(outliers)
            if length(outliers) > 1000, outliers = outliers(randperm(length(outliers), 1000)); end
            x_jit = center + (rand(size(outliers)) - 0.5) * 0.4;
            scatter(x_jit, outliers, 12, Color_Main, 'filled', 'MarkerFaceAlpha', 0.4);
        end
        % 箱体
        patch([center-0.25, center+0.25, center+0.25, center-0.25], [q1, q1, q3, q3], ...
              [0.96 0.96 0.96], 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
        plot([center-0.25, center+0.25], [med, med], 'Color', Color_Main, 'LineWidth', 2);
        plot([center, center], [p2_5, q1], 'k-', 'LineWidth', 1.2);
        plot([center, center], [q3, p97_5], 'k-', 'LineWidth', 1.2);
    end
end

% --- 优化三：均值收敛折线 ---
if Opt_MeanCurve
    % 画轨迹连线
    plot(1:8, means_val, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
    % 画均值点
    scatter(1:8, means_val, 60, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    % 高亮起点和终点
    scatter(1, means_val(1), 60, Color_Main, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

%% 4. 坐标轴及样式美化
set(ax, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Times New Roman');
xlim([0.3, 8.7]);
xticks(1:8);
xticklabels({'Term 1', 'Term 2', 'Term 3', 'Term 4', 'Term 5', 'Term 6', 'Term 7', 'Term 8'});
xlabel('Taylor Expansion Terms', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Trend Contribution (mm/yr)', 'FontSize', 14, 'FontWeight', 'bold');

% --- 优化四：SymLog 自定义坐标轴标签 ---
if Opt_SymLog
    % 还原 Y 轴真实的物理读数标签
    title('Convergence of Trend (SymLog Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 自定义要显示的真实物理刻度 (例如: 0.001, 0.01, 0.1, 1)
    real_ticks = [-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1];
    mapped_ticks = trans(real_ticks);
    
    yticks(mapped_ticks);
    yticklabels(arrayfun(@(x) num2str(x), real_ticks, 'UniformOutput', false));
    ylim([trans(-2), trans(2)]); % 限制显示范围
else
    title(sprintf('Convergence of %s (%s) Trend', dataset, state), 'FontSize', 16, 'FontWeight', 'bold');
    % 为了防止 1、2 阶太大导致高阶看不见，智能设置 YLim
    max_v = max(abs(prctile(Data_Terms{2}, [2.5, 97.5]))); 
    if max_v > 0, ylim([-max_v*1.5, max_v*1.5]); end
end

fprintf('>> 绘图完成！你可以修改脚本顶部的 4 个开关来查看不同效果。\n');
end