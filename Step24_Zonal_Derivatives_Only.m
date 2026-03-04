%% Step24_Zonal_Derivatives_Only.m
% 功能: 仅绘制【密度对温度的 N 阶导数】纬向平均剖面图
% 特点: 干净、光滑、中文标注、物理意义清晰
clear; clc; close all;

% ================= 1. 自动读取数据 =================
data_path = 'D:\work\EN4_analyses_c13_last20years\'; 
files = dir(fullfile(data_path, '**', '*.nc'));
if isempty(files), error('❌ 找不到数据文件'); end
file_path = fullfile(files(1).folder, files(1).name);
fprintf('>> 读取背景场: %s\n', files(1).name);

% 读取变量
lat = ncread(file_path, 'lat');
lon = ncread(file_path, 'lon');
depth = ncread(file_path, 'depth');
T_raw = ncread(file_path, 'temperature') - 273.15; 
S_raw = ncread(file_path, 'salinity');

% 经度对齐
[lon, sort_idx] = sort(mod(ncread(file_path,'lon'), 360));
T_raw = T_raw(sort_idx, :, :);
S_raw = S_raw(sort_idx, :, :);
lon(lon>180) = lon(lon>180)-360; 

% 网格生成
[LON, LAT, DEP] = ndgrid(lon, lat, depth);

% 检查 GSW
if isempty(which('gsw_rho'))
    error('❌ 请安装 GSW (TEOS-10) 工具箱');
end

% 计算压力 (若无函数则近似)
if ~isempty(which('gsw_p_from_z'))
    PRESS = gsw_p_from_z(-DEP, LAT);
else
    PRESS = abs(DEP);
end

% ================= 2. 核心：平滑与求导 =================
fprintf('>> 正在计算并平滑 (请稍候)...\n');

% 计算 SA, CT
SA = gsw_SA_from_SP(S_raw, PRESS, LON, LAT);
CT = gsw_CT_from_pt(SA, T_raw);

% 【关键】3D 高斯平滑 (去除网格噪音，保证画图光滑)
% 窗口 [7, 7, 3] 对应 [经度, 纬度, 深度]
CT_Smooth = smooth3(CT, 'gaussian', [7 7 3], 1.2);
SA_Smooth = smooth3(SA, 'gaussian', [7 7 3], 1.2);

% 定义求导函数
func_rho = @(ct_in) gsw_rho(SA_Smooth, ct_in, PRESS);
h = 0.2; % 差分步长

% 计算 1-3 阶导数
% 1阶: 热膨胀
D1 = (func_rho(CT_Smooth+h) - func_rho(CT_Smooth-h)) / (2*h);
% 2阶: 非线性/混合增密 (主导项)
D2 = (func_rho(CT_Smooth+h) - 2*func_rho(CT_Smooth) + func_rho(CT_Smooth-h)) / h^2;
% 3阶: 符号翻转项 (解释 T3 反号)
D3 = (func_rho(CT_Smooth+2*h) - 2*func_rho(CT_Smooth+h) + 2*func_rho(CT_Smooth-h) - func_rho(CT_Smooth-2*h)) / (2*h^3);

% ================= 3. 纬向平均 =================
% 沿经度(dim=1)平均
D1_Zonal = squeeze(nanmean(D1, 1)); 
D2_Zonal = squeeze(nanmean(D2, 1));
D3_Zonal = squeeze(nanmean(D3, 1));

% ================= 4. 绘图 (只画这一张) =================
figure('Position', [100, 100, 1500, 500], 'Color', 'w');
CN_Font = 'Microsoft YaHei'; 
[LAT_Grid, DEP_Grid] = meshgrid(lat, depth);

% 数据列表与标题
Data_List = {D1_Zonal, D2_Zonal, D3_Zonal};
Titles = {
    '(a) 一阶导数 \partial\rho/\partial T (热膨胀)', ...
    '(b) 二阶导数 \partial^2\rho/\partial T^2 (混合增密)', ...
    '(c) 三阶导数 \partial^3\rho/\partial T^3 (符号翻转)'
};

for i = 1:3
    subplot(1, 3, i);
    
    Data = Data_List{i};
    
    % 画填色图 (注意转置)
    pcolor(LAT_Grid, DEP_Grid, Data'); 
    shading interp; hold on;
    set(gca, 'YDir', 'reverse'); % 深度向下
    
    % 【关键】只画 0 等值线 (黑实线)
    contour(LAT_Grid, DEP_Grid, Data', [0 0], 'k', 'LineWidth', 2.0);
    
    % 设置色标 (红蓝平衡)
    % 如果有 cmocean 用 cmocean，没有则用自定义红蓝
    try
        colormap(gca, cmocean('balance'));
    catch
        % 自定义简易红蓝条
        R = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1)];
        B = [ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
        mycmap = [R; B];
        colormap(gca, mycmap);
    end
    
    % 自动调整色标范围 (对称)
    limit = prctile(abs(Data(:)), 98); 
    caxis([-limit, limit]);
    colorbar;
    
    % 标注与美化
    title(Titles{i}, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', CN_Font);
    xlabel('纬度 (Latitude)', 'FontSize', 12, 'FontName', CN_Font);
    ylabel('深度 (m)', 'FontSize', 12, 'FontName', CN_Font);
    ylim([0, 2000]); % 聚焦上层 2000m
    xlim([-80, 80]);
    
    % 在图上标注正负 (方便老师看)
    if i == 2
        text(0, 1000, '负值区域 (-)', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold', 'Horiz', 'center', 'FontName', CN_Font);
    elseif i == 3
        text(0, 500, '正值区域 (+)', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold', 'Horiz', 'center', 'FontName', CN_Font);
    end
end

sgtitle('海水密度对温度的 N 阶导数：纬向平均分布', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', CN_Font);
fprintf('>> 绘图完成！\n');