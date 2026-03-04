%% 趋势空间分布图 [StdRef Version]
%% Map_Trend_EN4.m
% 功能: EN4 数据趋势空间分布图 (显著性打点 + 全经纬度显示版)
% 修改点: 强制所有 8 个子图都显示经纬度刻度。
clear; clc; close all;

%% 1. 加载数据
files = {'EN4_TSLA_Terms_1to8_StdRef.mat'};
input_file = '';
DataDir = 'D:\work\EN4_TSLA_Terms';
for i = 1:length(files)
    if exist(files{i}, 'file')
        input_file = files{i};
    elseif exist(fullfile(DataDir, files{i}), 'file')
        input_file = fullfile(DataDir, files{i});
    end
    if ~isempty(input_file), break; end
end
if isempty(input_file)
    % 兼容旧文件名
    fallback_files = {'EN4_TSLA_Terms_1to15_StdRef.mat', 'EN4_TSLA_Terms_1to15_Safe.mat'};
    for i = 1:length(fallback_files)
        if exist(fallback_files{i}, 'file'), input_file = fallback_files{i}; break; end
        if exist(fullfile(DataDir, fallback_files{i}), 'file'), input_file = fullfile(DataDir, fallback_files{i}); break; end
    end
end
if isempty(input_file), error('❌ 找不到数据文件'); end
fprintf('>> 加载数据: %s ...\n', input_file);
load(input_file); 

% 经纬度修正
if ~exist('lon', 'var'), if exist('Lon','var'), lon=Lon; elseif exist('longitude','var'), lon=longitude; end; end
if ~exist('lat', 'var'), if exist('Lat','var'), lat=Lat; elseif exist('latitude','var'), lat=latitude; end; end
if max(lon) > 180, lon(lon > 180) = lon(lon > 180) - 360; end
[lon, sort_idx] = sort(lon);
TSLA_AllOrders = TSLA_AllOrders(sort_idx, :, :, :);
[nlon, nlat, ntime, max_order] = size(TSLA_AllOrders);
[Lon_Grid, Lat_Grid] = meshgrid(lon, lat);

% 计算 Trend 和 Significance
Map_Trend = nan(nlat, nlon, 8);
Map_Sig = nan(nlat, nlon, 8); 
t_vec = time_axis(:);
fprintf('>> 计算趋势及显著性检验 (p < 0.05)...\n');
for n = 1:8
    slice = squeeze(TSLA_AllOrders(:,:,:,n));
    temp_trend = nan(nlon, nlat);
    temp_sig = nan(nlon, nlat);
    
    for i=1:nlon
        for j=1:nlat
            ts = squeeze(slice(i,j,:));
            idx_valid = ~isnan(ts);
            if sum(idx_valid) < length(ts) * 0.7, continue; end 
            
            [p_coeff, S] = polyfit(t_vec(idx_valid), ts(idx_valid), 1);
            
            % StdRef 通常输出 m，所以 * 1000 转 mm
            temp_trend(i,j) = p_coeff(1) * 1000; 
            
            % t-test
            y_fit = polyval(p_coeff, t_vec(idx_valid));
            resid = ts(idx_valid) - y_fit;
            n_val = sum(idx_valid);
            X = [ones(n_val,1) t_vec(idx_valid)];
            XTX_inv = inv(X' * X);
            se_slope = sqrt(var(resid) * XTX_inv(2,2));
            t_stat = p_coeff(1) / se_slope;
            p_val = 2 * (1 - tcdf(abs(t_stat), n_val - 2));
            
            if p_val < 0.05
                temp_sig(i,j) = 1;
            end
        end
    end
    Map_Trend(:,:,n) = temp_trend';
    Map_Sig(:,:,n) = temp_sig';
end

% 智能修正：如果值过大，取消 * 1000
if max(abs(Map_Trend(:))) > 500
    Map_Trend = Map_Trend / 1000;
end

%% ==========================================
%% 绘图 (显著性叠加 + 全刻度)
%% ==========================================
figure('Position', [100, 50, 1200, 1400], 'Color', 'w'); 
% 稍微增大 TileSpacing 以防止刻度标签重叠
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% 配色定义
MyBlue = [0.05 0.25 0.6];   MyLightBlue = [0.6 0.8 1]; MyWhite = [1 1 1];
MyLightRed = [1 0.7 0.6];   MyRed = [0.7 0.05 0.05];
ColorLocs = [0, 0.25, 0.5, 0.75, 1]; 
ColorPoints = [MyBlue; MyLightBlue; MyWhite; MyLightRed; MyRed];
xx = linspace(0, 1, 256);
MyDivergentCmap = [interp1(ColorLocs, ColorPoints(:,1), xx, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,2), xx, 'pchip')', ...
                   interp1(ColorLocs, ColorPoints(:,3), xx, 'pchip')'];

for term_idx = 1:8
    nexttile;
    m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
    
    % 1. 绘制趋势
    data = Map_Trend(:,:,term_idx);
    m_pcolor(Lon_Grid, Lat_Grid, data); shading flat;
    hold on;
    
    % 2. 绘制显著性点
    sig_mask = Map_Sig(:,:,term_idx);
    skip = 3; 
    [row, col] = find(sig_mask(1:skip:end, 1:skip:end) == 1);
    if ~isempty(row)
        sig_lon = Lon_Grid(1:skip:end, 1:skip:end);
        sig_lat = Lat_Grid(1:skip:end, 1:skip:end);
        m_plot(sig_lon(sig_mask(1:skip:end, 1:skip:end)==1), ...
               sig_lat(sig_mask(1:skip:end, 1:skip:end)==1), ...
               'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
    end
    
    % 3. 海岸线
    m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
    
    % 4. 网格样式 (核心修改：移除所有隐藏逻辑，全部显示)
    % -----------------------------------------------------------
    m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
           'linewidth', 0.5, 'fontsize', 8); 
    % -----------------------------------------------------------
    
    colormap(gca, MyDivergentCmap);
    
    % 智能取整刻度
    raw_max = prctile(abs(data(:)), 99.8);
    if isempty(raw_max) || raw_max < 1e-12, raw_max = 1e-12; end
    exponent = floor(log10(raw_max)); 
    fraction = raw_max / 10^exponent; 
    if fraction <= 2, nice_base = 2; elseif fraction <= 5, nice_base = 5; else, nice_base = 10; end
    final_limit = nice_base * 10^exponent;
    caxis([-final_limit, final_limit]); 
    
    title(sprintf('Term %d', term_idx), 'FontSize', 12, 'FontWeight', 'bold');
    
    % 色棒
    cb = colorbar; 
    cb.FontSize = 8;
    title(cb, 'mm/yr', 'FontSize', 8); 
    cb.Ticks = linspace(-final_limit, final_limit, 5);
end

sgtitle('\fontname{Microsoft YaHei}EN4 (StdRef) 1-8 阶独立项趋势分布 (点区域 p < 0.05)', 'FontSize', 16, 'FontWeight', 'bold');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
fprintf('>> 绘图完成！每个子图均已显示经纬度刻度。\n');
