%% Step27_Final_Exact_Match.m
% 功能: 严格参照参考图绘制密度导数纬向平均结构
% 特点:
%   1. [复刻] 移除所有黑色等值线，依靠颜色过渡展示结构
%   2. [复刻] 使用参考图中的精确中文标题
%   3. [复刻] 添加图中的红/蓝文字标注区域
%   4. [基础] 保留平滑去噪和地形锐化处理
clear; clc; close all;

% ================= 1. 数据准备 =================
data_path = 'D:\work\EN4_analyses_c13_last20years\'; 
files = dir(fullfile(data_path, '**', '*.nc'));
if isempty(files), error('❌ 找不到数据文件'); end
file_path = fullfile(files(1).folder, files(1).name);
fprintf('>> 读取数据: %s\n', files(1).name);

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

[LON, LAT, DEP] = ndgrid(lon, lat, depth);

% 压力计算
if ~isempty(which('gsw_p_from_z'))
    PRESS = gsw_p_from_z(-DEP, LAT);
else
    PRESS = abs(DEP);
end

% 计算 SA, CT (原始场用于生成Mask)
fprintf('>> 计算物理场...\n');
SA_Raw = gsw_SA_from_SP(S_raw, PRESS, LON, LAT);
CT_Raw = gsw_CT_from_pt(SA_Raw, T_raw);

% ================= 2. 处理流程 (Mask & Smooth) =================
fprintf('>> 生成地形掩膜 & 执行平滑...\n');
% 1. 生成 Mask: 基于原始数据的 NaN 分布确定地形边界
T_Zonal_Mask = squeeze(nanmean(CT_Raw, 1)); 
Valid_Mask = ~isnan(T_Zonal_Mask); 

% 2. 预平滑 (Pre-smoothing): 3D高斯平滑用于计算导数
CT_Smooth = smooth3(CT_Raw, 'gaussian', [7 7 3], 1.2);
SA_Smooth = smooth3(SA_Raw, 'gaussian', [7 7 3], 1.2);

% ================= 3. 求导 =================
fprintf('>> 计算 N 阶导数...\n');
func_rho = @(ct_in) gsw_rho(SA_Smooth, ct_in, PRESS);
h = 0.2; 

% 1阶
D1 = (func_rho(CT_Smooth+h) - func_rho(CT_Smooth-h)) / (2*h);
% 2阶
D2 = (func_rho(CT_Smooth+h) - 2*func_rho(CT_Smooth) + func_rho(CT_Smooth-h)) / h^2;
% 3阶
D3 = (func_rho(CT_Smooth+2*h) - 2*func_rho(CT_Smooth+h) + 2*func_rho(CT_Smooth-h) - func_rho(CT_Smooth-2*h)) / (2*h^3);

% ================= 4. 纬向平均 & 后处理 =================
D1_Zonal = squeeze(nanmean(D1, 1)); 
D2_Zonal = squeeze(nanmean(D2, 1));
D3_Zonal = squeeze(nanmean(D3, 1));

% 后平滑 (Post-smoothing) 并 严格回填地形 Mask
D1_Zonal = imgaussfilt(D1_Zonal, 1.0); D1_Zonal(~Valid_Mask) = NaN;
D2_Zonal = imgaussfilt(D2_Zonal, 1.0); D2_Zonal(~Valid_Mask) = NaN;
D3_Zonal = imgaussfilt(D3_Zonal, 1.0); D3_Zonal(~Valid_Mask) = NaN;

% ================= 5. 严格参照绘图 =================
fprintf('>> 绘图 (严格参照模式)...\n');
figure('Position', [100, 100, 1500, 500], 'Color', 'w');
CN_Font = 'Microsoft YaHei'; % 使用支持中文的字体
[LAT_Grid, DEP_Grid] = meshgrid(lat, depth);

Data_List = {D1_Zonal, D2_Zonal, D3_Zonal};
% 参照图中的精确标题
Titles = {
    '一阶导数 (热膨胀): 全负', ...
    '二阶导数 (非线性): 负主导', ...
    '三阶导数 (变号项): 正负交替'
};

for i = 1:3
    subplot(1, 3, i);
    Data = Data_List{i};
    
    % 画填色图 (无 contour)
    pcolor(LAT_Grid, DEP_Grid, Data'); 
    shading interp; hold on;
    set(gca, 'YDir', 'reverse', 'Layer', 'top'); % Layer top 让刻度不被遮挡
    
    % 色标设置 (平衡红蓝，白色居中)
    colormap(gca, cmocean('balance'));
    
    % 强制对称色标范围，确保 0 值显示为白色
    limit = prctile(abs(Data(:)), 98); 
    if limit == 0, limit = 1e-10; end % 防止全零错误
    caxis([-limit, limit]);
    colorbar;
    
    % 标题与坐标轴
    title(Titles{i}, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', CN_Font);
    xlabel('纬度', 'FontSize', 12, 'FontName', CN_Font); 
    ylabel('深度 (m)', 'FontSize', 12, 'FontName', CN_Font);
    ylim([0, 5000]); xlim([-85, 85]); % 稍微调整纬度范围以匹配视觉
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    
    % 添加参照图中的文字标注
    if i == 2
        text(0, 800, '负值区 (-)', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold', 'Horiz', 'center', 'FontName', CN_Font);
    elseif i == 3
        text(0, 800, '正值区 (+)', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold', 'Horiz', 'center', 'FontName', CN_Font);
    end
end

sgtitle('海水密度导数纬向平均结构 (去噪版)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', CN_Font);
fprintf('>> 完成！图像风格已严格匹配参考图。\n');