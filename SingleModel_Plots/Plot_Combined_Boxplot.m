function Plot_Combined_Boxplot(state)
% 功能: 使用双/三 Y 轴的大一统溯源箱线图
% 设计: T1 使用左侧主 Y 轴；T2 及以后的高阶项根据 state 动态应用多轴叠加以防止刻度挤压。

clc; close all;
fprintf('========== 绘制统一溯源箱线图 ==========\n');

%% 1. 数据加载与准备
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
datasets = {'EN4', 'IAP', 'Ishii'};
% 顶刊溯源配色
colors = {
    [0.902, 0.294, 0.208], % EN4 (樱桃红)
    [0.302, 0.455, 0.690], % IAP (钴蓝)
    [0.000, 0.627, 0.408]  % Ishii (薄荷绿)
};

if strcmp(state, 'Average')
    max_terms = 3; fprintf('>> 检测到平均态，显示前 %d 阶...\n', max_terms);
else
    max_terms = 8; fprintf('>> 检测到标准态，显示全部 %d 阶...\n', max_terms);
end

Data_Trend = cell(3, max_terms);
for i = 1:length(datasets)
    trend_file = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', datasets{i}, state));
    if ~exist(trend_file, 'file'), error('❌ 找不到文件: %s', trend_file); end
    load(trend_file, 'trend_TSLA');
    for n = 1:max_terms
        slice = trend_TSLA(:,:,n);
        Data_Trend{i, n} = slice(~isnan(slice(:)));
    end
end

%% 2. 预设置参数与绘图布局
BoxWidth = 0.5; ScatterWidth = 0.4;
BoxBgColor = [0.92 0.92 0.92]; 

figure('Position', [100, 100, 1100, 550], 'Color', 'w', 'Name', 'Unified Boxplot');
MainPos = [0.08, 0.15, 0.82, 0.75];

if strcmp(state, 'Average')
    %% == 平均态：双 Y 轴分布 ==
    Pos_G = 1:max_terms;
    X_Limit = [0.4, max_terms + 0.6];
    for i = 1:max_terms, labels{i} = sprintf('T%d', i); end
    
    max_L = 0; max_R = 0;
    for n = 1:max_terms
        all_data = [Data_Trend{1, n}; Data_Trend{2, n}; Data_Trend{3, n}];
        if isempty(all_data), continue; end
        val = max(abs(prctile(all_data, [0.2, 99.8])));
        if n == 1, max_L = val; else, max_R = max(max_R, val); end
    end
    max_L = max_L * 1.2; if max_L == 0, max_L = 1; end
    max_R = max_R * 1.2; if max_R == 0, max_R = 1; end
    
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for term_idx = 1:max_terms
        center = Pos_G(term_idx);
        if term_idx == 1, yyaxis(ax1, 'left'); else yyaxis(ax1, 'right'); end
        draw_unified_box(ax1, center, Data_Trend(:, term_idx), colors, BoxWidth, ScatterWidth, BoxBgColor);
    end
    
    % 左轴修饰 (T1)
    yyaxis(ax1, 'left');
    ylim(ax1, [-max_L, max_L]);
    ax1.YColor = [0.15 0.15 0.15];
    ylabel(ax1, 'T1 Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14);
    yline(ax1, 0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    
    % 右轴修饰 (T2~T3)
    yyaxis(ax1, 'right');
    ylim(ax1, [-max_R, max_R]);
    ax1.YColor = [0.2 0.4 0.7];
    ylabel(ax1, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14, 'Color', [0.2 0.4 0.7]);
    
    set(ax1, 'XLim', X_Limit, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(ax1, 'XTick', Pos_G, 'XTickLabel', labels, 'FontWeight', 'bold');
    xline(ax1, 1.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    
    % == 图例 ==
    yyaxis(ax1, 'left'); 
    add_common_legend(ax1, BoxBgColor, colors, datasets);
    
    xlabel(ax1, 'Taylor Expansion Terms', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax1, sprintf('Unified Boxplot with Dual Y-axis (%%s)', state), 'FontSize', 16, 'FontWeight', 'bold');
    
else
    %% == 标准态：3 Y 轴联动分布 ==
    Pos_G1 = [1, 2, 3];
    Pos_G2 = [4.2, 5.2]; 
    Pos_G3 = [6.4, 7.4, 8.4];
    X_Limit = [0.4, 9.0];
    All_Ticks = [Pos_G1, Pos_G2, Pos_G3];
    All_Labels = {'T1','T2','T3', 'T4','T5', 'T6','T7','T8'};
    
    max_1 = 0; max_2 = 0; max_3 = 0;
    for n = 1:max_terms
        all_data = [Data_Trend{1, n}; Data_Trend{2, n}; Data_Trend{3, n}];
        if isempty(all_data), continue; end
        val = max(abs(prctile(all_data, [0.2, 99.8])));
        if ismember(n, 1:3), max_1 = max(max_1, val); end
        if ismember(n, 4:5), max_2 = max(max_2, val); end
        if ismember(n, 6:8), max_3 = max(max_3, val); end
    end
    max_1 = max_1 * 1.2; if max_1 == 0, max_1 = 1; end
    max_2 = max_2 * 1.2; if max_2 == 0, max_2 = 1; end
    max_3 = max_3 * 1.2; if max_3 == 0, max_3 = 1; end
    
    C_L1 = [0.15 0.15 0.15];
    C_L2 = [0.2 0.4 0.7]; 
    C_L3 = [0.000 0.627 0.408]; % 第三轴颜色（薄荷绿）
    
    % --- Layer 1 (T1-T3) ---
    ax1 = axes('Position', MainPos); hold(ax1, 'on');
    for k = 1:3
        draw_unified_box(ax1, Pos_G1(k), Data_Trend(:, k), colors, BoxWidth, ScatterWidth, BoxBgColor);
    end
    ylim(ax1, [-max_1, max_1]);
    set(ax1, 'XLim', X_Limit, 'YColor', C_L1, 'Box', 'off', 'XTick', [], 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax1, 'Dominant Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14, 'Color', C_L1);
    
    % --- Layer 2 (T4-T5) ---
    ax2 = axes('Position', MainPos); hold(ax2, 'on');
    for k = 1:2
        draw_unified_box(ax2, Pos_G2(k), Data_Trend(:, k+3), colors, BoxWidth, ScatterWidth, BoxBgColor);
    end
    ylim(ax2, [-max_2, max_2]);
    set(ax2, 'Color', 'none', 'XLim', X_Limit, 'YAxisLocation', 'right', 'YColor', C_L2, 'Box', 'off', 'XTick', [], 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax2, 'Secondary Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14, 'Color', C_L2);
    
    % --- Layer 3 (T6-T8) ---
    Offset_Ratio = 0.08; 
    Pos3 = MainPos; Pos3(3) = MainPos(3) * (1 + Offset_Ratio);
    ax3 = axes('Position', Pos3); hold(ax3, 'on');
    Scale_Factor = Pos3(3) / MainPos(3); 
    Real_X_Limit = X_Limit * Scale_Factor; 
    
    for k = 1:3
        draw_unified_box(ax3, Pos_G3(k), Data_Trend(:, k+5), colors, BoxWidth, ScatterWidth, BoxBgColor);
    end
    ylim(ax3, [-max_3, max_3]);
    set(ax3, 'Color', 'none', 'XLim', Real_X_Limit, 'YAxisLocation', 'right', 'YColor', C_L3, 'Box', 'off', 'XTick', [], 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel(ax3, 'High-Order Trend (mm/yr)', 'FontWeight', 'bold', 'FontSize', 14, 'Color', C_L3);
    
    % --- 分割线与基准线 ---
    yline(ax1, 0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5); 
    xline(ax1, 3.6, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    xline(ax1, 5.8, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    
    set(ax1, 'XTick', All_Ticks, 'XTickLabel', All_Labels, 'FontWeight', 'bold');
    xlabel(ax1, 'Taylor Expansion Terms', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax1, sprintf('Unified Boxplot with Triple Y-axis (%%s)', state), 'FontSize', 16, 'FontWeight', 'bold');
    
    % --- 同步联动 ---
    Ratio2 = max_2 / max_1; Ratio3 = max_3 / max_1;
    linkaxes([ax1, ax2], 'x'); 
    hZoom = zoom(gcf); hZoom.Motion = 'both'; hZoom.Enable = 'on';
    hZoom.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    hPan = pan(gcf); hPan.Motion = 'both';
    hPan.ActionPostCallback = @(s,e) SyncAllAxes(ax1, ax2, ax3, Ratio2, Ratio3, Scale_Factor);
    
    % == 图例 ==
    add_common_legend(ax1, BoxBgColor, colors, datasets);
end

fprintf('>> 多 Y 轴绘图完成！\n');
end

%% ============================================================
%% 内部函数
%% ============================================================
function draw_unified_box(ax, center, data_cells, colors, BoxWidth, ScatterWidth, BoxBgColor)
    all_data = [data_cells{1}; data_cells{2}; data_cells{3}];
    if isempty(all_data), return; end
    
    q1 = prctile(all_data, 25); q3 = prctile(all_data, 75);
    med_val = median(all_data); mean_val = mean(all_data);
    low_w = prctile(all_data, 2.5); up_w = prctile(all_data, 97.5);
    
    % 画离群点
    for i = 1:length(data_cells)
        curr_data = data_cells{i};
        idx_out = curr_data < low_w | curr_data > up_w;
        outliers = curr_data(idx_out);
        
        if ~isempty(outliers)
            if length(outliers) > 10000, outliers = outliers(randperm(length(outliers), 10000)); end
            % 按照参考图样式：均匀分布
            x_jit = center + (rand(size(outliers)) - 0.5) * ScatterWidth;
            
            % 按照参考图样式：增加黑色边缘，取消透明，尺寸适中
            scatter(ax, x_jit, outliers, 15, colors{i}, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.9);
        end
    end
    
    % 画须线
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w], 'k-', 'LineWidth', 1.2);
    cap_w = BoxWidth * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.5);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w], 'k-', 'LineWidth', 1.5);
    
    % 画箱体
    patch(ax, [center-BoxWidth/2, center+BoxWidth/2, center+BoxWidth/2, center-BoxWidth/2], ...
        [q1, q1, q3, q3], BoxBgColor, 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot(ax, [center-BoxWidth/2, center+BoxWidth/2], [med_val, med_val], 'k-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.2);
end

function add_common_legend(ax, BoxBgColor, colors, datasets)
    h_box  = patch(ax, NaN, NaN, BoxBgColor, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.85);
    h_med  = plot(ax, NaN, NaN, 'k-', 'LineWidth', 2);
    h_mean = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.2);
    h_out1 = scatter(ax, NaN, NaN, 50, colors{1}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.9);
    h_out2 = scatter(ax, NaN, NaN, 50, colors{2}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.9);
    h_out3 = scatter(ax, NaN, NaN, 50, colors{3}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.9);

    legend(ax, [h_box, h_med, h_mean, h_out1, h_out2, h_out3], ...
           {'IQR Box (25%-75%)', 'Median', 'Mean', [datasets{1} ' Outliers'], [datasets{2} ' Outliers'], [datasets{3} ' Outliers']}, ...
           'Location', 'northeast', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'none', 'FontName', 'Times New Roman');
end

function SyncAllAxes(ax1, ax2, ax3, R2, R3, Factor)
    try
        L = max(abs(ax1.YLim));
        ax1.YLim = [-L, L];
        ax2.YLim = [-L*R2, L*R2];
        ax3.YLim = [-L*R3, L*R3];
        ax3.XLim = ax1.XLim * Factor;
    catch
    end
end