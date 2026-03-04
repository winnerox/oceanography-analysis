%% Compare_Avg_vs_Std_Fixed.m
% =========================================================================
% 功能：对比 [平均态] vs [标准态] 的全球时间序列与差异场
% 【核心升级】：加入了二阶温盐交叉项 (Cross Term)
%
% 逻辑：
%   1. Avg State = Sum(Terms 1-3) + Cross_Avg
%   2. Std State = [Sum(Terms 1-8) - Mean] + Cross_Std_TSLA
% =========================================================================
clear; clc; close all;

%% 1. 文件路径配置
% -------------------------------------------------------------------------
File_Avg   = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms1to8_Average.mat';
File_Std   = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms_1to8_StdRef.mat';
File_Cross = 'D:\work\EN4_TSLA_Terms\EN4_CrossTerm_DualMode.mat'; % 新增

if ~exist(File_Avg, 'file') || ~exist(File_Std, 'file') || ~exist(File_Cross, 'file')
    error('❌ 找不到输入文件，请检查路径配置！');
end

%% 2. 加载数据
% -------------------------------------------------------------------------
fprintf('>> [1/6] 加载数据 (含交叉项)...\n');
D_Avg   = load(File_Avg);
D_Std   = load(File_Std);
D_Cross = load(File_Cross);

% 检查维度
fprintf('   Main Data Size: %s\n', mat2str(size(D_Avg.TSLA_AllOrders)));

%% 3. 网格对齐
% -------------------------------------------------------------------------
fprintf('>> [2/6] 网格对齐...\n');
% A. 主项对齐 (4D)
[Avg_Data, Avg_Lon, Avg_Lat] = align_grid(D_Avg.TSLA_AllOrders, D_Avg.lon, D_Avg.lat);
[Std_Data, Std_Lon, Std_Lat] = align_grid(D_Std.TSLA_AllOrders, D_Std.lon, D_Std.lat);

% B. 交叉项对齐 (3D)
[Cross_Avg_Data, ~, ~] = align_grid(D_Cross.Cross_Avg, D_Cross.lon, D_Cross.lat);
[Cross_Std_Data, ~, ~] = align_grid(D_Cross.Cross_Std_TSLA, D_Cross.lon, D_Cross.lat);

% C. 维度一致性检查
[nlon, nlat, ntime, ~] = size(Avg_Data);
if size(Cross_Avg_Data, 3) ~= ntime
    warning('⚠️ 时间维度不一致！主文件:%d, 交叉项:%d。正在截断...', ntime, size(Cross_Avg_Data,3));
    min_t = min(ntime, size(Cross_Avg_Data, 3));
    Avg_Data = Avg_Data(:,:,1:min_t,:);
    Std_Data = Std_Data(:,:,1:min_t,:);
    Cross_Avg_Data = Cross_Avg_Data(:,:,1:min_t);
    Cross_Std_Data = Cross_Std_Data(:,:,1:min_t);
    D_Avg.time_axis = D_Avg.time_axis(1:min_t);
end

%% 4. 执行求和与组合 (Main + Cross)
% -------------------------------------------------------------------------
fprintf('>> [3/6] 组合最终 TSLA (Main + Cross)...\n');

% --- A. 平均态组合 ---
% 1. 主项求和 (1-3阶)
Sum_Avg_Main = nansum(Avg_Data(:,:,:,1:3), 4);
Mask_Avg = all(isnan(Avg_Data(:,:,:,1:3)), 4);
Sum_Avg_Main(Mask_Avg) = NaN;

% 2. 加上交叉项
Final_Avg = Sum_Avg_Main + Cross_Avg_Data;


% --- B. 标准态组合 ---
% 1. 主项求和 (1-8阶)
Sum_Std_Total = nansum(Std_Data(:,:,:,1:8), 4);
Mask_Std = all(isnan(Std_Data(:,:,:,1:8)), 4);
Sum_Std_Total(Mask_Std) = NaN;

% 2. 去平均 (Convert to Anomaly)
Mean_Field_Main = nanmean(Sum_Std_Total, 3);
Sum_Std_Main_Anomaly = Sum_Std_Total - Mean_Field_Main;

% 3. 加上交叉项 (Cross_Std_Data 已经在生成时去过平均了)
Final_Std = Sum_Std_Main_Anomaly + Cross_Std_Data;

% 输出均值检查
fprintf('   Final Avg Mean: %.4f mm\n', nanmean(Final_Avg(:)));
fprintf('   Final Std Mean: %.4f mm\n', nanmean(Final_Std(:)));

%% 5. 计算差值 (Diff)
% -------------------------------------------------------------------------
fprintf('>> [4/6] 计算差异场 (Std_Final - Avg_Final)...\n');

Diff_Map_3D = Final_Std - Final_Avg;

% A. Bias
Mean_Diff_Map = nanmean(Diff_Map_3D, 3);
% B. RMSE
RMS_Diff_Map = sqrt(nanmean(Diff_Map_3D.^2, 3));

fprintf('   差异范围: [%.2f, %.2f] mm\n', min(Diff_Map_3D(:)), max(Diff_Map_3D(:)));

%% 6. 计算全球平均时间序列
% -------------------------------------------------------------------------
fprintf('>> [5/6] 计算全球加权平均曲线...\n');
[W_Lon, W_Lat] = meshgrid(1:length(Avg_Lon), cosd(Avg_Lat));
Weight_Grid = W_Lat'; 

calc_weighted_mean = @(data_3d, w_2d) ...
    squeeze(nansum(nansum(data_3d .* w_2d, 1), 2) ./ ...
            nansum(nansum(~isnan(data_3d) .* w_2d, 1), 2));

Global_TS_Avg  = calc_weighted_mean(Final_Avg, Weight_Grid);
Global_TS_Std  = calc_weighted_mean(Final_Std, Weight_Grid);
Global_TS_Diff = Global_TS_Std - Global_TS_Avg;

t_vec = D_Avg.time_axis;

%% 7. 绘图：空间分布差异
% -------------------------------------------------------------------------
fprintf('>> [6/6] 开始绘图...\n');
[Lon_Grid, Lat_Grid] = meshgrid(Avg_Lon, Avg_Lat);
Cmap_Diff = m_colmap('diverging'); % 使用 M_Map 自带颜色

figure('Position', [100, 100, 1000, 800], 'Color', 'w');
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Bias ---
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, Mean_Diff_Map'); shading flat;
m_coast('patch', [.9 .9 .9], 'edgecolor', [.5 .5 .5]);
m_grid('linestyle', ':', 'tickdir', 'out', 'linewidth', 0.5);
colormap(gca, Cmap_Diff); 
caxis([-5, 5]); 
cb = colorbar; title(cb, 'mm');
title({'\fontname{Microsoft YaHei}TSLA 系统差异 (Bias)', ...
       '\fontsize{10}\color{gray}[Std(1-8)+Cross] - [Avg(1-3)+Cross]'}, 'FontSize', 12, 'FontWeight', 'bold');

% --- RMSE ---
nexttile;
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Grid, Lat_Grid, RMS_Diff_Map'); shading flat;
m_coast('patch', [.9 .9 .9], 'edgecolor', [.5 .5 .5]);
m_grid('linestyle', ':', 'tickdir', 'out', 'linewidth', 0.5);
colormap(gca, 'parula'); 
caxis([0, 5]); 
cb = colorbar; title(cb, 'mm');
title({'\fontname{Microsoft YaHei}波动差异幅度 (RMSE)', ...
       '\fontsize{10}\color{gray}包含交叉项后的总误差'}, 'FontSize', 12, 'FontWeight', 'bold');

sgtitle('\fontname{Microsoft YaHei}全物理量对比: 平均态 vs 标准态', 'FontSize', 16, 'FontWeight', 'bold');

%% 8. 绘图：全球时间序列对比
% -------------------------------------------------------------------------
figure('Position', [150, 150, 1000, 400], 'Color', 'w');
plot(t_vec, Global_TS_Std, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Std State (1-8 + Cross)');
hold on;
plot(t_vec, Global_TS_Avg, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Avg State (1-3 + Cross)');
plot(t_vec, Global_TS_Diff, 'k-', 'LineWidth', 1, 'DisplayName', 'Difference');
yline(0, 'k:', 'LineWidth', 0.8);

legend('Location', 'best', 'FontSize', 10);
xlabel('Year'); ylabel('Global Mean TSLA (mm)');
title('\fontname{Microsoft YaHei}全球 TSLA 时间序列终极对比 (含二阶交叉项)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

fprintf('>> ✅ 全流程对比完成！\n');

%% ============================================================
%% 内部函数
%% ============================================================
function [data_out, lon_out, lat_out] = align_grid(data_in, lon_in, lat_in)
    % 兼容 3D 和 4D 数据的网格对齐
    if max(lon_in) > 180, lon_in(lon_in > 180) = lon_in(lon_in > 180) - 360; end
    [lon_out, sort_idx] = sort(lon_in);
    lat_out = lat_in;
    
    if ndims(data_in) == 4
        data_out = data_in(sort_idx, :, :, :);
    else
        data_out = data_in(sort_idx, :, :);
    end
end