%% Analysis_CrossTerms.m
% =========================================================================
% 功能：混合项 (Cross Terms) 量级对比分析
% 对比：Average State vs Standard Reference State
% 
% 输入：
%   1. EN4_Terms_S_Cross_Avg_1to3.mat (Avg State Cross: ST, SST, STT)
%   2. EN4_CrossTerms_Order8_Mismatch.mat (StdRef Cross: HO_Cross_Std_TSLA)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\CrossTerms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> 加载数据...\n');

% Avg State (分项: ST, SST, STT)
f_avg = fullfile(DataDir, 'EN4_Terms_S_Cross_Avg_1to3.mat');
if ~exist(f_avg, 'file'), error('❌ 请先运行 Calc_Cross_and_S_Terms_Avg.m!'); end
D_Avg = load(f_avg, 'Cross_ST', 'Cross_SST', 'Cross_STT', 'lon', 'lat');

% StdRef State (总和: HO_Cross_Std_TSLA)
f_std = fullfile(DataDir, 'EN4_CrossTerms_Order8_Mismatch.mat');
if ~exist(f_std, 'file')
    f_std = fullfile(DataDir, 'EN4_CrossTerms_HighOrder_DualMode.mat');
end
if exist(f_std, 'file')
    D_Std = load(f_std);
    % 实际变量名是 HO_Cross_Std_TSLA (Total, not per-term)
    if isfield(D_Std, 'HO_Cross_Std_TSLA')
        Cross_Std_Total = D_Std.HO_Cross_Std_TSLA;
    else
        Cross_Std_Total = nan;
    end
    if isfield(D_Std, 'Lon'), lon_std = D_Std.Lon; else, lon_std = D_Avg.lon; end
    if isfield(D_Std, 'Lat'), lat_std = D_Std.Lat; else, lat_std = D_Avg.lat; end
else
    error('❌ 找不到标准态混合项数据文件.');
end

%% 3. 计算统计量 (RMS Magnitude)
fprintf('>> 计算量级对比...\n');

% 计算全球 RMS 时间平均量级的函数
calc_mag = @(x) nanmean(sqrt(nanmean(x.^2, 3)), 'all'); 

% Avg: 分项
Mag_ST_Avg = calc_mag(D_Avg.Cross_ST);
Mag_SST_Avg = calc_mag(D_Avg.Cross_SST);
Mag_STT_Avg = calc_mag(D_Avg.Cross_STT);
Mag_Total_Avg = calc_mag(D_Avg.Cross_ST + D_Avg.Cross_SST + D_Avg.Cross_STT);

% Std: 只有总和
Mag_Total_Std = calc_mag(Cross_Std_Total);

fprintf('\n=== 混合项量级对比 (RMS mm) ===\n');
fprintf('项         | 平均态   | 标准态\n');
fprintf('----------------------------------\n');
fprintf('ST (2阶)   | %.4f   | (含于总和)\n', Mag_ST_Avg);
fprintf('SST (3阶)  | %.4f   | (含于总和)\n', Mag_SST_Avg);
fprintf('STT (3阶)  | %.4f   | (含于总和)\n', Mag_STT_Avg);
fprintf('----------------------------------\n');
fprintf('总混合项    | %.4f   | %.4f\n', Mag_Total_Avg, Mag_Total_Std);
fprintf('----------------------------------\n');

%% 4. 绘图 (柱状图 - 总量对比)
figure('Position', [100, 100, 500, 400], 'Color', 'w');
bar([Mag_Total_Avg, Mag_Total_Std]);
xticklabels({'平均态', '标准态'});
ylabel('全球 RMS 量级 (mm)');
title('混合项总量级对比');
grid on;
saveas(gcf, fullfile(OutputDir, 'CrossTerm_Comparison_Bar.png'));

%% 5. 空间分布对比 (地图)
% 选取 总混合项 (Total Cross) 进行对比
figure('Position', [100, 100, 1000, 500], 'Color', 'w');
t = tiledlayout(1, 2, 'Padding', 'compact');

% Calc Trend
Cross_Avg_Total = D_Avg.Cross_ST + D_Avg.Cross_SST + D_Avg.Cross_STT;
Trend_Avg = calc_trend(Cross_Avg_Total);
Trend_Std = calc_trend(Cross_Std_Total);

% Avg
nexttile; 
map_plot(D_Avg.lon, D_Avg.lat, Trend_Avg', '平均态: 混合项趋势', [-0.2 0.2]);

% Std
nexttile;
map_plot(lon_std, lat_std, Trend_Std', '标准态: 混合项趋势', [-0.2 0.2]);

sgtitle('混合项 (所有阶数之和) 趋势对比');
saveas(gcf, fullfile(OutputDir, 'CrossTerm_Total_Map.png'));

fprintf('>> 分析完成！图像已保存至 %s\n', OutputDir);

%% Helper
function t = calc_trend(data)
    [Nx, Ny, Nt] = size(data);
    t = nan(Nx, Ny);
    x = (1:Nt)';
    for i=1:Nx
        for j=1:Ny
            y = squeeze(data(i,j,:));
            if sum(~isnan(y)) > Nt*0.8
                p = polyfit(x, y, 1);
                t(i,j) = p(1) * 12; % mm/yr (monthly data)
            end
        end
    end
end

function map_plot(lon, lat, data, tit, clim)
    % 输入: lon [Nx], lat [Ny], data [Ny x Nx] (已转置)
    lon = lon(:)'; lat = lat(:)';
    
    % 经度转换：如果是 0-360 格式，转为 -180 to 180
    if max(lon) > 180
        shift_idx = lon > 180;
        lon(shift_idx) = lon(shift_idx) - 360;
        [lon, sort_idx] = sort(lon);
        data = data(:, sort_idx); % 重排列 (lon 是第二维)
    end
    
    % 创建网格
    [LON, LAT] = meshgrid(lon, lat);
    
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    m_pcolor(LON, LAT, data); shading flat;
    m_coast('patch', [.9 .9 .9], 'edgecolor', 'k');
    m_grid('linestyle', 'none', 'tickdir', 'out');
    colorbar; caxis(clim);
    title(tit);
    colormap(m_colmap('diverging'));
end

