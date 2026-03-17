function Plot_Exact_SSLA_Difference()
% 功能: 绘制 EN4, IAP, Ishii 三套数据集下，Average态与StdRef态精确解(SSLA)、热膨胀(TSLA)、盐收缩(HSLA)的差异分布图
% 特性: 3x3 全景并排对比，独立网格适配，揭示计算底层一致性及机制拆解偏移

clc; close all;
fprintf('========== 🌊 启动参考态一致性与拆解差异 (SSLA/TSLA/HSLA) 验证引擎 ==========\n');

%% 1. 初始化设置
datasets = {'EN4', 'IAP', 'Ishii'};
trend_dir = 'D:\work\Task_Convergence\Trend_Results';

% 存储插值/处理后差异及独立网格的 Cell 数组
diff_ssla_maps = cell(1, 3);
diff_tsla_maps = cell(1, 3);
diff_hsla_maps = cell(1, 3);
Lon_Grids = cell(1, 3);
Lat_Grids = cell(1, 3);

% 存储各行全局最大的残差，用于统一色棒
global_max_ssla = 0;
global_max_tsla = 0;
global_max_hsla = 0;

%% 2. 批量处理三套数据
for d = 1:length(datasets)
    ds = datasets{d};
    fprintf('>> 正在处理 [%s] 数据集...\n', ds);
    
    avg_file = fullfile(trend_dir, sprintf('%s_Average_Trends.mat', ds));
    std_file = fullfile(trend_dir, sprintf('%s_StdRef_Trends.mat', ds));
    
    if ~exist(avg_file, 'file') || ~exist(std_file, 'file')
        error('❌ 找不到 %s 的趋势数据文件，请检查路径！', ds);
    end
    
    d_avg = load(avg_file);
    d_std = load(std_file);
    
    % 提取精确解
    ssla_avg = d_avg.trend_SSLA(:,:,1);
    ssla_std = d_std.trend_SSLA(:,:,1);
    
    % 提取并聚合 TSLA 和 HSLA (1-8阶求和)
    tsla_avg = sum(d_avg.trend_TSLA, 3, 'omitnan');
    tsla_std = sum(d_std.trend_TSLA, 3, 'omitnan');
    hsla_avg = sum(d_avg.trend_HSLA, 3, 'omitnan');
    hsla_std = sum(d_std.trend_HSLA, 3, 'omitnan');
    
    % 恢复陆地掩膜 (以精确解的掩膜为基准)
    mask = ~isnan(ssla_avg);
    ssla_std(~mask) = NaN;
    tsla_avg(~mask) = NaN; tsla_std(~mask) = NaN;
    hsla_avg(~mask) = NaN; hsla_std(~mask) = NaN;
    
    % 智能单位换算 (m -> mm)
    if max(abs(ssla_avg(:))) < 0.2
        ssla_avg = ssla_avg * 1000; ssla_std = ssla_std * 1000;
        tsla_avg = tsla_avg * 1000; tsla_std = tsla_std * 1000;
        hsla_avg = hsla_avg * 1000; hsla_std = hsla_std * 1000;
    end
    
    % 计算各分量的参考态差异 (Average - StdRef)
    diff_ssla = ssla_avg - ssla_std;
    diff_tsla = tsla_avg - tsla_std;
    diff_hsla = hsla_avg - hsla_std;
    
    % --- 经纬度对齐处理 ---
    Lon = double(d_avg.Lon(:))';
    Lat = double(d_avg.Lat(:))';
    nLon = length(Lon);
    
    % 转置安检强制 [Lat x Lon]
    if size(diff_ssla, 1) == nLon
        diff_ssla = diff_ssla';
        diff_tsla = diff_tsla';
        diff_hsla = diff_hsla';
    end
    
    % 缝合线对齐 (0~360 -> -180~180)
    if max(Lon) > 180
        Lon(Lon > 180) = Lon(Lon > 180) - 360;
        [Lon, sort_idx] = sort(Lon);
        diff_ssla = diff_ssla(:, sort_idx);
        diff_tsla = diff_tsla(:, sort_idx);
        diff_hsla = diff_hsla(:, sort_idx);
    end
    
    % 保存独立专属经纬度网格与数据
    [Lon_Grids{d}, Lat_Grids{d}] = meshgrid(Lon, Lat);
    diff_ssla_maps{d} = diff_ssla;
    diff_tsla_maps{d} = diff_tsla;
    diff_hsla_maps{d} = diff_hsla;
    
    % 寻找 99.8% 极值，避免孤立极值噪点拉毁各行色棒
    global_max_ssla = max(global_max_ssla, prctile(abs(diff_ssla(:)), 99.8));
    global_max_tsla = max(global_max_tsla, prctile(abs(diff_tsla(:)), 99.8));
    global_max_hsla = max(global_max_hsla, prctile(abs(diff_hsla(:)), 99.8));
end

%% 3. 计算极值与色棒范围限制
% SSLA (精确解极值应无限趋近于 0)
if global_max_ssla < 1e-15, global_max_ssla = 1e-15; end
exp_ssla = floor(log10(global_max_ssla)); frac_ssla = global_max_ssla / 10^exp_ssla;
if frac_ssla <= 2, nb = 2; elseif frac_ssla <= 5, nb = 5; else, nb = 10; end
limit_ssla = nb * 10^exp_ssla;

% TSLA (热膨胀差异通常较大)
if global_max_tsla < 1e-12, global_max_tsla = 1e-12; end
exp_tsla = floor(log10(global_max_tsla)); frac_tsla = global_max_tsla / 10^exp_tsla;
if frac_tsla <= 2, nb = 2; elseif frac_tsla <= 5, nb = 5; else, nb = 10; end
limit_tsla = nb * 10^exp_tsla;

% HSLA (盐收缩差异)
if global_max_hsla < 1e-12, global_max_hsla = 1e-12; end
exp_hsla = floor(log10(global_max_hsla)); frac_hsla = global_max_hsla / 10^exp_hsla;
if frac_hsla <= 2, nb = 2; elseif frac_hsla <= 5, nb = 5; else, nb = 10; end
limit_hsla = nb * 10^exp_hsla;

%% 4. 绘制 3x3 顶刊级验证长卷
fprintf('\n>> 开始渲染 3x3 空间机制拆解阵列...\n');
figure('Position', [50, 50, 1800, 1300], 'Color', 'w', 'Name', 'Components State Diff Maps');
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

% =================== 第 1 行: Exact SSLA ===================
for d = 1:3
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grids{d}, Lat_Grids{d}, diff_ssla_maps{d}); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_ssla, limit_ssla]);
    
    title(sprintf('%s Exact SSLA Diff (Avg-Std)', datasets{d}), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
    cb = colorbar('Location', 'southoutside'); title(cb, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
    cb.Ticks = linspace(-limit_ssla, limit_ssla, 5);
end

% =================== 第 2 行: TSLA ===================
for d = 1:3
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grids{d}, Lat_Grids{d}, diff_tsla_maps{d}); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_tsla, limit_tsla]);
    
    title(sprintf('%s TSLA Diff (Avg-Std)', datasets{d}), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.8 0.2 0.2]);
    cb = colorbar('Location', 'southoutside'); title(cb, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
    cb.Ticks = linspace(-limit_tsla, limit_tsla, 5);
end

% =================== 第 3 行: HSLA ===================
for d = 1:3
    nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(Lon_Grids{d}, Lat_Grids{d}, diff_hsla_maps{d}); shading flat; hold on;
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
    colormap(gca, MyDivergentCmap); caxis([-limit_hsla, limit_hsla]);
    
    title(sprintf('%s HSLA Diff (Avg-Std)', datasets{d}), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.1 0.3 0.7]);
    cb = colorbar('Location', 'southoutside'); title(cb, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
    cb.Ticks = linspace(-limit_hsla, limit_hsla, 5);
end

% --- 总标题 ---
sgtitle('Difference of SSLA, TSLA, HSLA (Average - StdRef)', 'FontSize', 22, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 🎉 3x3 全景验证图绘制完成！\n');

%% 5. 保存结果
save_dir = 'D:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fig_name = 'All_Components_State_Difference.jpg';
filename = fullfile(save_dir, fig_name);
saveas(gcf, filename, 'jpg'); 
fprintf('>> 已保存: %s\n', filename);

end