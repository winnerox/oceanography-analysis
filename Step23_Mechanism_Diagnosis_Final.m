%% Step50_HighOrder_8_Levels.m
% 功能: 利用高阶引擎计算并绘制 1 到 8 阶密度导数 (d^n rho / dT^n)
% 核心: 循环调用 Engine -> 智能布局 -> 深度探索状态方程非线性
clear; clc; close all;

% ================= 1. 数据准备 (3D 平均场) =================
data_path = 'D:\work\EN4_analyses_c13_last20years\';
files = dir(fullfile(data_path, '**', '*.nc'));
if isempty(files), error('❌ 找不到数据文件'); end

num_read = min(length(files)); 
fprintf('>> 读取 %d 个文件计算 3D 平均场...\n', num_read);

f1 = fullfile(files(1).folder, files(1).name);
lat = ncread(f1, 'lat'); lon = ncread(f1, 'lon'); depth = ncread(f1, 'depth');
[lon, sort_idx] = sort(mod(lon, 360)); lon(lon>180) = lon(lon>180)-360;

sz = [length(lon), length(lat), length(depth)];
Sum_T = zeros(sz); Count_T = zeros(sz);
Sum_S = zeros(sz); Count_S = zeros(sz);

for i = 1:num_read
    fname = fullfile(files(i).folder, files(i).name);
    t_raw = ncread(fname, 'temperature') - 273.15;
    s_raw = ncread(fname, 'salinity');
    t_raw = t_raw(sort_idx, :, :); s_raw = s_raw(sort_idx, :, :);
    
    valid_t = ~isnan(t_raw); t_raw(~valid_t) = 0;
    Sum_T = Sum_T + t_raw; Count_T = Count_T + valid_t;
    valid_s = ~isnan(s_raw); s_raw(~valid_s) = 0;
    Sum_S = Sum_S + s_raw; Count_S = Count_S + valid_s;
end
T_mean = Sum_T ./ Count_T; T_mean(Count_T == 0) = NaN;
S_mean = Sum_S ./ Count_S; S_mean(Count_S == 0) = NaN;

% ================= 2. 批量计算 1-8 阶导数 =================
fprintf('>> 准备计算 1 到 8 阶导数...\n');
[LON, LAT, DEP] = ndgrid(lon, lat, depth);
if ~isempty(which('gsw_p_from_z')), PRESS = gsw_p_from_z(-DEP, LAT); else, PRESS = abs(DEP); end

% 计算基础场
SA = gsw_SA_from_SP(S_mean, PRESS, LON, LAT);
CT = gsw_CT_from_pt(SA, T_mean);

% 平滑 (对高阶求导至关重要，否则高阶全是噪点)
% 随着阶数增加，建议稍微加大一点平滑力度，或者保持适中
SA = smooth3(SA, 'gaussian', [5 5 3], 0.8);
CT = smooth3(CT, 'gaussian', [5 5 3], 0.8);

% 初始化引擎
try
    eng = TEOS10_HighOrder_Engine();
catch
    error('请确保 TEOS10_HighOrder_Engine.m 在当前路径下');
end

% 存储结果的元胞数组
Results = cell(1, 8);
Max_Order = 8;

fprintf('   正在计算:\n');
for n = 1:Max_Order
    fprintf('   -> 第 %d 阶导数 (d^%d rho / dT^%d)... ', n, n, n);
    t_start = tic;
    
    % 调用引擎: calculate_mixed(SA, CT, p, n_T, n_S)
    % 这里 n_T = n, n_S = 0 (纯温度导数)
    rho_deriv = eng.calculate_mixed(SA, CT, PRESS, n, 0);
    
    % 立即计算纬向平均以节省内存
    Results{n} = squeeze(nanmean(rho_deriv, 1));
    
    fprintf('耗时 %.2f 秒\n', toc(t_start));
end

% ================= 3. 绘图 (2x4 布局) =================
fprintf('>> 绘图...\n');
figure('Position', [50, 50, 1800, 800], 'Color', 'w');
t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

mycolor = cmocean('balance');

for n = 1:Max_Order
    nexttile;
    
    data = Results{n};
    
    % 绘图 (Teacher Style)
    h_img = imagesc(lat, depth, data');
    set(gca, 'Color', [0.8 0.8 0.8]); % 浅灰背景
    set(h_img, 'AlphaData', ~isnan(data')); % 透明遮罩
    set(gca, 'YDir', 'reverse');
    
    colormap(gca, mycolor);
    
    % --- 智能色标 ---
    valid_data = data(~isnan(data));
    low = prctile(valid_data, 2);
    high = prctile(valid_data, 98);
    
    if n <= 2
        % 1阶和2阶通常全是负的，使用非对称色标
        % 稍微扩宽一点范围让颜色柔和
        span = high - low;
        caxis([low - span*0.1, high + span*0.2]);
    else
        % 3阶及以上通常是震荡的 (正负交替)，强制对称色标
        limit = max(abs([low, high]));
        caxis([-limit, limit] * 1.1); % 乘 1.1 柔化
    end
    % ----------------
    
    % 装饰
    title(sprintf('Order %d (\\partial^%d\\rho / \\partial T^%d)', n, n, n), 'FontSize', 12, 'FontWeight', 'bold');
    if n > 4, xlabel('Latitude'); end
    if mod(n, 4) == 1, ylabel('Depth (m)'); end
    
    cb = colorbar;
    cb.Label.String = sprintf('kg m^{-3} °C^{-%d}', n);
end

title(t, '密度对温度的 1 到 8 阶导数纬向分布 (TEOS-10 Engine)', 'FontSize', 16, 'FontWeight', 'bold');