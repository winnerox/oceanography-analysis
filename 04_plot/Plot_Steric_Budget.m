function Plot_Steric_Budget(state)
% 功能: 绘制全球比容海平面趋势收支闭合图 (Steric Budget)
% 包含: 分组对比柱状图 + 批量生成 4 张 1x5空间机制拆解与残差验证图 (EN4/IAP/Ishii/Ensemble)
% 修复: 独立尺寸安检，彻底解决矩阵经纬度错位导致的 interp2 崩溃问题
% 参数: state = 'StdRef' 或 'Average'

clc; close all;
fprintf('========== 🌊 启动终极收支机制拆解 (Budget Breakdown) 分析 [%s] ==========\n', state);

%% 1. 数据准备与全局均值计算
datasets = {'EN4', 'IAP', 'Ishii'};
trend_dir = 'D:\work\Task_Convergence\Trend_Results';

gm_SSLA = zeros(1, 4); % 1-3为单源, 4为Ensemble
gm_TSLA = zeros(1, 4);
gm_HSLA = zeros(1, 4);
gm_Cross = zeros(1, 4);

% 用于收集 Ensemble 的三维矩阵
target_lon = 0.5:1:359.5; target_lat = -89.5:1:89.5;
[Target_Lon_Grid, Target_Lat_Grid] = meshgrid(target_lon, target_lat);

Ens_SSLA = []; Ens_TSLA = []; Ens_HSLA = []; Ens_Cross = []; Ens_Res = [];

for d = 1:3
    ds = datasets{d};
    f_trend = fullfile(trend_dir, sprintf('%s_%s_Trends.mat', ds, state));
    f_cross = fullfile(trend_dir, sprintf('%s_%s_Cross_Trends.mat', ds, state));
    
    if ~exist(f_trend, 'file') || ~exist(f_cross, 'file')
        error('❌ 找不到数据文件，请确保已经跑完了 Trend 和 Cross 的计算脚本。');
    end
    
    t_main = load(f_trend);
    t_cross = load(f_cross);
    
    % 提取并求和各分量 (泰勒阶数求和)
    T_SSLA  = t_main.trend_SSLA(:,:,1); % 精确解
    T_TSLA  = sum(t_main.trend_TSLA, 3, 'omitnan');
    T_HSLA  = sum(t_main.trend_HSLA, 3, 'omitnan');
    T_Cross = sum(t_cross.trend_Cross, 3, 'omitnan');
    
    % 获取经纬度基准
    curr_lon = double(t_main.Lon(:))'; 
    curr_lat = double(t_main.Lat(:))';
    nLon = length(curr_lon); 
    nLat = length(curr_lat);
    
    % ==============================================================
    % 🚨 核心修复：分别对每个矩阵进行独立安检！
    % 确保所有矩阵绝对被强制规范为 [Lat x Lon] 格式，杜绝插值和掩膜崩溃
    % ==============================================================
    if size(T_SSLA, 1) == nLon, T_SSLA = T_SSLA'; end
    if size(T_TSLA, 1) == nLon, T_TSLA = T_TSLA'; end
    if size(T_HSLA, 1) == nLon, T_HSLA = T_HSLA'; end
    if size(T_Cross, 1) == nLon, T_Cross = T_Cross'; end
    
    % 恢复陆地 NaN 掩膜 (现在尺寸绝对对齐，安全赋值)
    mask = ~isnan(T_SSLA);
    T_TSLA(~mask) = NaN; T_HSLA(~mask) = NaN; T_Cross(~mask) = NaN;
    
    % 单位换算 (m -> mm)
    if max(abs(T_SSLA(:))) < 0.2
        T_SSLA = T_SSLA * 1000; T_TSLA = T_TSLA * 1000;
        T_HSLA = T_HSLA * 1000; T_Cross = T_Cross * 1000;
    end
    
    % --- 🚀 计算全球面积加权平均 ---
    [~, Curr_Lat_Grid] = meshgrid(curr_lon, curr_lat);
    weights = cosd(Curr_Lat_Grid); % 纬度余弦加权
    weights(~mask) = NaN;
    
    gm_SSLA(d) = nansum(T_SSLA(mask) .* weights(mask)) / nansum(weights(mask));
    gm_TSLA(d) = nansum(T_TSLA(mask) .* weights(mask)) / nansum(weights(mask));
    gm_HSLA(d) = nansum(T_HSLA(mask) .* weights(mask)) / nansum(weights(mask));
    gm_Cross(d) = nansum(T_Cross(mask) .* weights(mask)) / nansum(weights(mask));
    
    fprintf('>> [%s] 全球均值:\n   Exact = %.5f | TSLA = %.5f | HSLA = %.5f | Cross = %.5f\n', ...
        ds, gm_SSLA(d), gm_TSLA(d), gm_HSLA(d), gm_Cross(d));
        
    % --- 空间插值对齐，准备堆叠 3D 矩阵 ---
    if min(curr_lon) < 0
        curr_lon(curr_lon < 0) = curr_lon(curr_lon < 0) + 360;
        [curr_lon, sort_idx] = sort(curr_lon);
        T_SSLA = T_SSLA(:, sort_idx); T_TSLA = T_TSLA(:, sort_idx);
        T_HSLA = T_HSLA(:, sort_idx); T_Cross = T_Cross(:, sort_idx);
    end
    curr_lon = [curr_lon(end)-360, curr_lon, curr_lon(1)+360];
    T_SSLA = cat(2, T_SSLA(:,end), T_SSLA, T_SSLA(:,1));
    T_TSLA = cat(2, T_TSLA(:,end), T_TSLA, T_TSLA(:,1));
    T_HSLA = cat(2, T_HSLA(:,end), T_HSLA, T_HSLA(:,1));
    T_Cross = cat(2, T_Cross(:,end), T_Cross, T_Cross(:,1));
    
    [Curr_Lon_Grid, Curr_Lat_Grid] = meshgrid(curr_lon, curr_lat);
    SSLA_interp = interp2(Curr_Lon_Grid, Curr_Lat_Grid, T_SSLA, Target_Lon_Grid, Target_Lat_Grid, 'linear');
    TSLA_interp = interp2(Curr_Lon_Grid, Curr_Lat_Grid, T_TSLA, Target_Lon_Grid, Target_Lat_Grid, 'linear');
    HSLA_interp = interp2(Curr_Lon_Grid, Curr_Lat_Grid, T_HSLA, Target_Lon_Grid, Target_Lat_Grid, 'linear');
    Cross_interp = interp2(Curr_Lon_Grid, Curr_Lat_Grid, T_Cross, Target_Lon_Grid, Target_Lat_Grid, 'linear');
    
    % 核心机制验证：精确解减去所有偏导和的残差
    Res_interp = SSLA_interp - (TSLA_interp + HSLA_interp + Cross_interp);
    
    Ens_SSLA = cat(3, Ens_SSLA, SSLA_interp);
    Ens_TSLA = cat(3, Ens_TSLA, TSLA_interp);
    Ens_HSLA = cat(3, Ens_HSLA, HSLA_interp);
    Ens_Cross = cat(3, Ens_Cross, Cross_interp);
    Ens_Res = cat(3, Ens_Res, Res_interp);
end

% 🚀 融合生成 Ensemble 集合平均全球均值
gm_SSLA(4) = mean(gm_SSLA(1:3));
gm_TSLA(4) = mean(gm_TSLA(1:3));
gm_HSLA(4) = mean(gm_HSLA(1:3));
gm_Cross(4) = mean(gm_Cross(1:3));

%% 2. 绘制图表 1: 全球平均收支对比柱状图 (分组并列 Grouped Bar)
fprintf('\n>> 绘制全球收支分组柱状图...\n');
figure('Position', [100, 100, 1100, 600], 'Color', 'w', 'Name', 'Global Budget Bar Chart');
ax = axes('Position', [0.08, 0.15, 0.88, 0.75]); hold on;

C_Exact = [0.4 0.4 0.4];
C_TSLA  = [0.902, 0.294, 0.208];
C_HSLA  = [0.302, 0.455, 0.690];
C_Cross = [0.850, 0.500, 0.100];

% 组合数据 4x4 (行: 数据集, 列: 分量)
y_data = [gm_SSLA', gm_TSLA', gm_HSLA', gm_Cross'];

% 绘制簇状柱形图
b = bar(y_data, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
b(1).FaceColor = C_Exact;
b(2).FaceColor = C_TSLA;
b(3).FaceColor = C_HSLA;
b(4).FaceColor = C_Cross;

% 添加完美的残差文本提示
for i = 1:4
    res = gm_SSLA(i) - (gm_TSLA(i) + gm_HSLA(i) + gm_Cross(i));
    y_pos = max([gm_SSLA(i), gm_TSLA(i), gm_HSLA(i)]) + 0.15;
    if y_pos < 0.2, y_pos = 0.2; end
    
    text(i, y_pos, sprintf('Closure Residual:\n%.1e', res), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.6 0.1 0.1]);
end

% 装饰与坐标轴
set(ax, 'XTick', 1:4, 'XTickLabel', {'EN4', 'IAP', 'Ishii', 'Ensemble Mean'}, ...
    'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'LineWidth', 1.2);
ylabel('Global Mean Trend (mm/yr)', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Global Mean Steric Sea Level Budget Breakdown (%s)', state), 'FontSize', 18, 'FontWeight', 'bold');
grid on; set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5);

legend({'Exact SSLA', 'TSLA (Total Thermal)', 'HSLA (Total Haline)', 'Cross (Total Nonlinear)'}, ...
    'Location', 'northwest', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'k');

%% 3. 批量绘制图表 2: 空间物理机制拆解图 (1x5 闭合验证长卷 - 针对每一套数据)
fprintf('>> 批量绘制空间机制闭合验证长卷 (EN4, IAP, Ishii, Ensemble)...\n');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

% 提前将矩阵经度对齐到 -180~180
Lon_E = target_lon; Lat_E = target_lat;
if max(Lon_E) > 180
    Lon_E(Lon_E > 180) = Lon_E(Lon_E > 180) - 360;
    [Lon_E, sort_idx] = sort(Lon_E);
    Ens_SSLA = Ens_SSLA(:, sort_idx, :);
    Ens_TSLA = Ens_TSLA(:, sort_idx, :);
    Ens_HSLA = Ens_HSLA(:, sort_idx, :);
    Ens_Cross = Ens_Cross(:, sort_idx, :);
    Ens_Res = Ens_Res(:, sort_idx, :);
end
[Lon_Grid_E, Lat_Grid_E] = meshgrid(Lon_E, Lat_E);

dataset_labels = {'EN4', 'IAP', 'Ishii', 'Ensemble Mean'};

% 🚀 循环绘制 4 张 1x5 空间图！
for k = 1:4
    figure('Position', [50, 100 + (k-1)*30, 2200, 450], 'Color', 'w', 'Name', sprintf('Spatial Breakdown Map %s', dataset_labels{k}));
    t = tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 提取当前数据集的切片
    if k == 4
        curr_SSLA = mean(Ens_SSLA, 3, 'omitnan');
        curr_TSLA = mean(Ens_TSLA, 3, 'omitnan');
        curr_HSLA = mean(Ens_HSLA, 3, 'omitnan');
        curr_Cross = mean(Ens_Cross, 3, 'omitnan');
        curr_Res = mean(Ens_Res, 3, 'omitnan');
    else
        curr_SSLA = Ens_SSLA(:,:,k);
        curr_TSLA = Ens_TSLA(:,:,k);
        curr_HSLA = Ens_HSLA(:,:,k);
        curr_Cross = Ens_Cross(:,:,k);
        curr_Res = Ens_Res(:,:,k);
    end
    
    % 核心色棒自适应
    max_main = max([prctile(abs(curr_SSLA(:)), 99.5), prctile(abs(curr_TSLA(:)), 99.5), prctile(abs(curr_HSLA(:)), 99.5)]);
    limit_main = ceil(max_main);
    if limit_main < 1, limit_main = 1; end

    % --------- (1) Exact SSLA (精确总解) ---------
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid_E, Lat_Grid_E, curr_SSLA); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_main, limit_main]);
    title('1. Exact SSLA', 'FontSize', 14, 'FontWeight', 'bold');
    cb1 = colorbar('Location', 'southoutside'); title(cb1, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');

    % --------- (2) TSLA (纯热膨胀) ---------
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid_E, Lat_Grid_E, curr_TSLA); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_main, limit_main]);
    title('2. TSLA (\Sigma 1-8)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', C_TSLA);
    cb2 = colorbar('Location', 'southoutside'); title(cb2, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');

    % --------- (3) HSLA (纯盐收缩) ---------
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid_E, Lat_Grid_E, curr_HSLA); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_main, limit_main]);
    title('3. HSLA (\Sigma 1-8)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.1 0.3 0.6]);
    cb3 = colorbar('Location', 'southoutside'); title(cb3, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');

    % --------- (4) Cross Terms (温盐混合非线性项) ---------
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid_E, Lat_Grid_E, curr_Cross); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); 

    if strcmp(state, 'Average'), pct = 99.0; else, pct = 99.5; end
    raw_max = prctile(abs(curr_Cross(:)), pct); 
    if isempty(raw_max) || isnan(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); fraction = raw_max / 10^exponent; 
    if fraction <= 2, nb = 2; elseif fraction <= 5, nb = 5; else, nb = 10; end
    limit_cross = nb * 10^exponent; caxis([-limit_cross, limit_cross]);

    title('4. Cross (\Sigma 1-28)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', C_Cross);
    cb4 = colorbar('Location', 'southoutside'); title(cb4, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
    cb4.Ticks = linspace(-limit_cross, limit_cross, 5);

    % --------- (5) Residual (严丝合缝的闭合证明) ---------
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grid_E, Lat_Grid_E, curr_Res); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); 

    % 残差极限锁定
    res_max = max(1e-12, prctile(abs(curr_Res(:)), 99.8)); 
    exponent_res = floor(log10(res_max)); 
    fraction_res = res_max / 10^exponent_res; 
    if fraction_res <= 2, nb_r = 2; elseif fraction_res <= 5, nb_r = 5; else, nb_r = 10; end
    limit_res = nb_r * 10^exponent_res; caxis([-limit_res, limit_res]);

    title('5. Residual (1 - 2 - 3 - 4)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.8 0.1 0.1]);
    cb5 = colorbar('Location', 'southoutside'); title(cb5, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
    cb5.Ticks = linspace(-limit_res, limit_res, 5);

    % 总标题明确标示数据来源
    sgtitle(sprintf('[%s] Spatial Budget Breakdown (%s)', dataset_labels{k}, state), 'FontSize', 18, 'FontWeight', 'bold');
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
end

fprintf('>> 🎉 完美拆解图全部绘制完成！\n');

%% 4. 保存结果
save_dir = 'D:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i); if isvalid(fig)
        fig_name = strrep(get(fig, 'Name'), ' ', '_');
        filename = fullfile(save_dir, sprintf('Budget_%s_%s.jpg', state, fig_name));
        saveas(fig, filename, 'jpg'); 
    end
end
end