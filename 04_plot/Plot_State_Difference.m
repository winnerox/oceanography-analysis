function Plot_State_Difference(dataset_name, target_comp)
% 功能: 绘制 Average 态与 StdRef 态 1-8 阶累加和的空间差异分布图
% 特性: 1x3 对比长卷 (Average Sum | StdRef Sum | Difference)
% 输入参数:
%   dataset_name: 数据集名称 (例如 'EN4', 'IAP', 'Ishii')
%   target_comp: 目标变量 ('TSLA', 'HSLA' 或 'Taylor_SSLA')，默认 'TSLA'

clc; close all;
fprintf('========== 🌊 启动参考态差异 (Average vs StdRef) 分析引擎 ==========\n');

if nargin < 2
    target_comp = 'TSLA'; % 默认对比热膨胀
end

%% 1. 加载两个状态的结果文件
trend_dir = 'D:\work\Task_Convergence\Trend_Results';
avg_file = fullfile(trend_dir, sprintf('%s_Average_Trends.mat', dataset_name));
std_file = fullfile(trend_dir, sprintf('%s_StdRef_Trends.mat', dataset_name));

if ~exist(avg_file, 'file'), error('❌ 找不到 Average 趋势文件: %s', avg_file); end
if ~exist(std_file, 'file'), error('❌ 找不到 StdRef 趋势文件: %s', std_file); end

fprintf('>> 正在加载 [%s] Average 态数据...\n', dataset_name);
d_avg = load(avg_file);
fprintf('>> 正在加载 [%s] StdRef 态数据...\n', dataset_name);
d_std = load(std_file);

%% 2. 提取 1-8 阶并求和
fprintf('>> 🎯 当前对比变量: 【%s】\n', target_comp);

if strcmpi(target_comp, 'TSLA')
    var_name = 'trend_TSLA';
elseif strcmpi(target_comp, 'HSLA')
    var_name = 'trend_HSLA';
elseif strcmpi(target_comp, 'Taylor_SSLA')
    % 如果是 SSLA，则合成 TSLA + HSLA
    var_name = 'Taylor_SSLA';
else
    error('❌ 不支持的变量: %s', target_comp);
end

% 提取并求和 (Average)
if strcmp(var_name, 'Taylor_SSLA')
    sum_avg = sum(d_avg.trend_TSLA, 3, 'omitnan') + sum(d_avg.trend_HSLA, 3, 'omitnan');
else
    sum_avg = sum(d_avg.(var_name), 3, 'omitnan');
end

% 提取并求和 (StdRef)
if strcmp(var_name, 'Taylor_SSLA')
    sum_std = sum(d_std.trend_TSLA, 3, 'omitnan') + sum(d_std.trend_HSLA, 3, 'omitnan');
else
    sum_std = sum(d_std.(var_name), 3, 'omitnan');
end

% 获取原始掩膜以恢复陆地 NaN
mask = ~isnan(d_avg.trend_SSLA(:,:,1));
sum_avg(~mask) = NaN;
sum_std(~mask) = NaN;

% 智能单位换算 (m -> mm)
if max(abs(sum_avg(:))) < 0.2
    sum_avg = sum_avg * 1000;
    sum_std = sum_std * 1000;
end

% 计算核心差异矩阵: Average - StdRef
diff_map = sum_avg - sum_std;

%% 3. 经纬度对齐处理
Lon = double(d_avg.Lon(:))';
Lat = double(d_avg.Lat(:))';
nLon = length(Lon);

% 转置安检
if size(sum_avg, 1) == nLon
    sum_avg = sum_avg';
    sum_std = sum_std';
    diff_map = diff_map';
end

% 缝合线对齐 (0~360 -> -180~180)
if max(Lon) > 180
    Lon(Lon > 180) = Lon(Lon > 180) - 360;
    [Lon, sort_idx] = sort(Lon);
    sum_avg = sum_avg(:, sort_idx);
    sum_std = sum_std(:, sort_idx);
    diff_map = diff_map(:, sort_idx);
end
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

%% 4. 绘制 1x3 顶刊级对比长卷
fprintf('>> 开始渲染 1x3 空间对比长卷...\n');
figure('Position', [50, 200, 1800, 500], 'Color', 'w', 'Name', sprintf('%s State Difference', target_comp));
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 顶刊发散型色盘
RdBu_ColorPoints = [0.1 0.2 0.5; 0.6 0.8 0.9; 1.0 1.0 1.0; 1.0 0.7 0.7; 0.6 0.1 0.1];
RdBu_Locs = [0, 0.3, 0.5, 0.7, 1]; xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(RdBu_Locs, RdBu_ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(RdBu_Locs, RdBu_ColorPoints(:,3), xx, 'pchip')'];

% 统一前两张图(主信号)的色棒极值
max_main = max(prctile(abs(sum_avg(:)), 99.5), prctile(abs(sum_std(:)), 99.5));
limit_main = ceil(max_main);
if limit_main < 1, limit_main = 1; end

% --------- (1) Average State Sum ---------
nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, sum_avg); shading flat; hold on;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
colormap(gca, MyDivergentCmap); caxis([-limit_main, limit_main]);
title('1. Average State (\Sigma 1-8 Terms)', 'FontSize', 14, 'FontWeight', 'bold');
cb1 = colorbar('Location', 'southoutside'); title(cb1, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');

% --------- (2) StdRef State Sum ---------
nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, sum_std); shading flat; hold on;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
colormap(gca, MyDivergentCmap); caxis([-limit_main, limit_main]);
title('2. StdRef State (\Sigma 1-8 Terms)', 'FontSize', 14, 'FontWeight', 'bold');
cb2 = colorbar('Location', 'southoutside'); title(cb2, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');

% --------- (3) Difference (Average - StdRef) ---------
nexttile; m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, diff_map); shading flat; hold on;
m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], 'linewidth', 0.5, 'fontsize', 8); 
colormap(gca, MyDivergentCmap); 

% 差异图拥有独立的、更精细的色棒限制，以放大细微的偏移模式
diff_max = prctile(abs(diff_map(:)), 99.5);
if isempty(diff_max) || diff_max < 1e-12, diff_max = 1e-12; end
exponent = floor(log10(diff_max)); fraction = diff_max / 10^exponent; 
if fraction <= 2, nb = 2; elseif fraction <= 5, nb = 5; else, nb = 10; end
limit_diff = nb * 10^exponent; 
caxis([-limit_diff, limit_diff]);

title('3. Difference (Average - StdRef)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1]);
cb3 = colorbar('Location', 'southoutside'); title(cb3, 'mm/yr', 'FontSize', 10, 'FontWeight', 'bold');
cb3.Ticks = linspace(-limit_diff, limit_diff, 5);

% --- 总标题 ---
sgtitle(sprintf('[%s] Effect of Reference State on %s Cumulative Sum', dataset_name, strtok(target_comp, '_')), ...
        'FontSize', 18, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('>> 🎉 差异对比图绘制完成！\n');

%% 5. 保存结果
save_dir = 'D:\work\Figures'; if ~exist(save_dir, 'dir'), mkdir(save_dir); end
figs = findall(0, 'Type', 'figure');
for i = 1:length(figs)
    fig = figs(i); if isvalid(fig)
        fig_name = strrep(get(fig, 'Name'), ' ', '_');
        filename = fullfile(save_dir, sprintf('Diff_%s_%s_%s.jpg', dataset_name, target_comp, fig_name));
        saveas(fig, filename, 'jpg'); 
    end
end
end