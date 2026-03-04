%% Spatial_Diff_Trend_Amp_RealSalt.m
% =========================================================================
% 功能：验证 [真实盐度背景下的 8 阶展开] 是否能让趋势图收敛
% 
% 核心逻辑：
%   Avg State: 依然用原来的 (1-3阶 + 混合项) -> 代表传统做法
%   Std State: 换成新的 (真实盐度背景 1-8阶) -> 代表你的改进做法
% 
% 预期结果：
%   深红色大斑块消失，全图变白（或仅剩细微的蓝色噪点）。
% =========================================================================
clear; clc; close all;

%% 1. 文件路径配置
% -------------------------------------------------------------------------
% 1. 平均态文件 (保持不变)
File_Avg   = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms1to8_Average.mat';
% 2. 平均态混合项 (保持不变)
File_Cross = 'D:\work\EN4_TSLA_Terms\EN4_CrossTerms_Order8_Mismatch.mat'; 

% 3. 【核心替换】标准态文件 -> 换成你刚刚跑完的“真实盐度版”
File_Std_New = 'D:\work\EN4_TSLA_Terms\EN4_TSLA_Terms_1to8_StdRef_RealSalt.mat';

if ~exist(File_Avg, 'file') || ~exist(File_Std_New, 'file')
    error('❌ 找不到文件，请检查路径！');
end

%% 2. 加载数据
fprintf('>> [1/5] 加载数据...\n');
D_Avg   = load(File_Avg);
D_Cross = load(File_Cross);
D_Std   = load(File_Std_New); % 加载新文件

% 网格对齐
[Avg_Main, Avg_Lon, Avg_Lat] = align_grid(D_Avg.TSLA_AllOrders, D_Avg.lon, D_Avg.lat);
[Cross_A, ~, ~]              = align_grid(D_Cross.HO_Cross_Avg, D_Cross.Lon, D_Cross.Lat);

% 【注意】新文件的 Std 数据
[Std_Main, Std_Lon, Std_Lat] = align_grid(D_Std.TSLA_AllOrders, D_Std.lon, D_Std.lat);

% 时间对齐
t_vec = D_Avg.time_axis(:);
min_t = min([size(Avg_Main,3), size(Std_Main,3), size(Cross_A,3)]);
Avg_Main = Avg_Main(:,:,1:min_t,:); 
Std_Main = Std_Main(:,:,1:min_t,:);
Cross_A  = Cross_A(:,:,1:min_t);
t_vec = t_vec(1:min_t);

[nlon, nlat, ntime, ~] = size(Avg_Main);

%% 3. 构造对比序列 (Total Budget)
fprintf('>> [2/5] 构造总量...\n');

% === A. 平均态组合 (Avg State) ===
% 逻辑：主项(1-3) + 混合项(2-3)
% 这是目前公认比较准的近似方法
TSLA_Avg_Total = nansum(Avg_Main(:,:,:,1:3), 4) + Cross_A;
TSLA_Avg_Total(all(isnan(Avg_Main(:,:,:,1:3)), 4)) = NaN;

% === B. 新标准态组合 (Real Salt Std State) ===
% 逻辑：直接把 1-8 阶加起来即可。
% 【关键点】：因为我们用了真实盐度，Delta_S = 0。
% 所以：不需要加 Cross Terms (混合项)，也不需要加 Salinity Terms (纯盐项)。
% 它本身就是包含所有非线性的热贡献。
TSLA_Std_Raw = nansum(Std_Main(:,:,:,1:8), 4); 

% 去平均 (Anomaly)
% 这一步至关重要，因为两个方法的零点定义不同，必须看 Anomaly
TSLA_Std_Total = TSLA_Std_Raw - nanmean(TSLA_Std_Raw, 3);
TSLA_Avg_Total = TSLA_Avg_Total - nanmean(TSLA_Avg_Total, 3);

%% 4. 计算趋势 (Trend)
fprintf('>> [3/5] 计算趋势差异...\n');
Trend_Avg = nan(nlat, nlon); 
Trend_Std = nan(nlat, nlon); 

for i = 1:nlon
    for j = 1:nlat
        % Avg
        y1 = squeeze(TSLA_Avg_Total(i,j,:));
        if sum(~isnan(y1)) > ntime*0.7
            p1 = polyfit(t_vec(~isnan(y1)), y1(~isnan(y1)), 1);
            Trend_Avg(j,i) = p1(1);
        end
        % Std
        y2 = squeeze(TSLA_Std_Total(i,j,:));
        if sum(~isnan(y2)) > ntime*0.7
            p2 = polyfit(t_vec(~isnan(y2)), y2(~isnan(y2)), 1);
            Trend_Std(j,i) = p2(1);
        end
    end
end

%% 5. 计算差异
Diff_Trend = Trend_Std - Trend_Avg;
fprintf('   Trend 差异范围: [%.4f, %.4f] mm/yr\n', min(Diff_Trend(:)), max(Diff_Trend(:)));

%% 6. 绘图
fprintf('>> [5/5] 绘图...\n');
[Lon_Plot, Lat_Plot] = meshgrid(Avg_Lon, Avg_Lat);

figure('Position', [100, 100, 1000, 600], 'Color', 'w');
m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
m_pcolor(Lon_Plot, Lat_Plot, Diff_Trend); shading flat;
m_coast('patch', [.9 .9 .9], 'edgecolor', [.5 .5 .5]); m_grid('linestyle',':');

colormap(m_colmap('diverging')); 
caxis([-0.5, 0.5]); % 保持这个严格的色标
cb = colorbar; title(cb, 'mm/yr');

title({'\fontname{Microsoft YaHei}最终闭合验证：真实盐度修正版', ...
       '\color{blue}New Std (Real Salt, Order 1-8) \color{black}vs \color{red}Avg State (Order 3)', ...
       '\fontsize{12}预期：深红色斑块应消失，整体趋于白色'}, ...
       'FontSize', 14, 'FontWeight', 'bold');

if max(abs(Diff_Trend(:))) < 0.2
    fprintf('>> 🎉 完美！差异非常小，证明两种方法殊途同归。\n');
else
    fprintf('>> 差异减小了，但如果还有残留，可能是“盐容海平面”的微小贡献差异。\n');
end

%% 辅助函数
function [data_out, lon_out, lat_out] = align_grid(data_in, lon_in, lat_in)
    if max(lon_in) > 180, lon_in(lon_in > 180) = lon_in(lon_in > 180) - 360; end
    [lon_out, sort_idx] = sort(lon_in);
    lat_out = lat_in;
    if ndims(data_in) == 4
        data_out = data_in(sort_idx, :, :, :);
    elseif ndims(data_in) == 3
        data_out = data_in(sort_idx, :, :);
    else
        data_out = data_in(sort_idx, :); 
    end
end