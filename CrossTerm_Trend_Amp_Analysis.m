%% Step3_HighOrder_Extended_FinalFixed.m
% =========================================================================
% 功能：计算【标准态】下的 2, 3, 4, 5 阶温盐混合偏导项影响
% 状态：【终极修复版】
%   1. 修复 shading flat 在匿名函数中报错的问题（改为本地函数）。
%   2. 修复 parfor 导致的编码报错（改为普通 for）。
%   3. 包含完整的 2-5 阶泰勒展开系数。
% =========================================================================
clear; clc; close all;

%% 1. 配置区域
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms_HighOrder';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
Years = 2005:2024; 
Months = 1:12;
MaxDepth = 2000;     
rho0 = 1035.0; 

% === 标准态定义 ===
Std_SA_Val = 35.0;  
Std_CT_Val = 0.0;   

if isempty(which('gsw_rho_high_order'))
    error('❌ 请确保 gsw_rho_high_order.m (或 fast版本) 在路径中！');
end

%% 2. 初始化网格
fprintf('>> [1/7] 初始化网格...\n');
SampleFile = '';
for y = Years
    FList = dir(fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y), '*.nc'));
    if ~isempty(FList), SampleFile = fullfile(FList(1).folder, FList(1).name); break; end
end
if isempty(SampleFile), error('❌ 找不到数据文件'); end

Lon_Orig = ncread(SampleFile, 'lon');
Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon_Orig, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);

Nz = length(CalcDepth);
dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

%% 3. 设置参考场
SA_Ref = Std_SA_Val * ones(size(P_3D));
CT_Ref = Std_CT_Val * ones(size(P_3D));

%% 4. 批量计算高阶导数系数
fprintf('>> [2/7] 计算 2-5 阶混合偏导数系数...\n');
[Nx, Ny, ~] = size(P_3D);
get_rho_deriv = @(m, n) reshape(gsw_rho_high_order.derivative(SA_Ref, CT_Ref, P_3D, m, n, 'S', 'T'), Nx, Ny, Nz);

% --- 2阶 ---
R_ST = get_rho_deriv(1,1);
% --- 3阶 ---
R_SST = get_rho_deriv(2,1); R_STT = get_rho_deriv(1,2);
% --- 4阶 ---
R_SSST = get_rho_deriv(3,1); R_SSTT = get_rho_deriv(2,2); R_STTT = get_rho_deriv(1,3);
% --- 5阶 ---
R_SSSST = get_rho_deriv(4,1); R_SSSTT = get_rho_deriv(3,2); R_SSTTT = get_rho_deriv(2,3); R_STTTT = get_rho_deriv(1,4);

%% 5. 遍历时间序列
fprintf('>> [3/7] 计算混合项贡献 (泰勒展开)... \n');
TotalFiles = length(Years)*12;

% Fig1 Data
TS_2_ST = zeros(Nx, Ny, TotalFiles, 'single');
TS_3_SST = zeros(Nx, Ny, TotalFiles, 'single'); TS_3_STT = zeros(Nx, Ny, TotalFiles, 'single');
% Fig2 Data (4阶)
TS_4_SSST = zeros(Nx, Ny, TotalFiles, 'single'); TS_4_SSTT = zeros(Nx, Ny, TotalFiles, 'single'); TS_4_STTT = zeros(Nx, Ny, TotalFiles, 'single');
% Fig2 Data (5阶)
TS_5_SSSST = zeros(Nx, Ny, TotalFiles, 'single'); TS_5_SSSTT = zeros(Nx, Ny, TotalFiles, 'single'); TS_5_SSTTT = zeros(Nx, Ny, TotalFiles, 'single'); TS_5_STTTT = zeros(Nx, Ny, TotalFiles, 'single');

time_vec = zeros(TotalFiles, 1);
idx = 0;

for y = Years
    for m = Months
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        
        Files = dir(fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y), sprintf('*%d%02d*.nc', y, m)));
        if isempty(Files), continue; end
        
        try
            T_Inst = double(squeeze(ncread(fullfile(Files(1).folder, Files(1).name), 'temperature', [1 1 1 1], [Inf Inf Nz 1])));
            S_Inst = double(squeeze(ncread(fullfile(Files(1).folder, Files(1).name), 'salinity', [1 1 1 1], [Inf Inf Nz 1])));
        catch
            continue; 
        end
        if nanmean(T_Inst(:)) > 100, T_Inst = T_Inst - 273.15; end
        
        SA_Inst = gsw_SA_from_SP(S_Inst, P_3D, LON_3D, LAT_3D);
        CT_Inst = gsw_CT_from_pt(SA_Inst, T_Inst);
        dS = SA_Inst - SA_Ref; 
        dT = CT_Inst - CT_Ref; 
        
        % 泰勒项计算 (带阶乘系数)
        t_st = R_ST .* dS .* dT;
        
        t_sst = (1/2) * R_SST .* (dS.^2) .* dT;
        t_stt = (1/2) * R_STT .* dS .* (dT.^2);
        
        t_ssst = (1/6) * R_SSST .* (dS.^3) .* dT;
        t_sstt = (1/4) * R_SSTT .* (dS.^2) .* (dT.^2);
        t_sttt = (1/6) * R_STTT .* dS .* (dT.^3);
        
        t_sssst = (1/24) * R_SSSST .* (dS.^4) .* dT;
        t_ssstt = (1/12) * R_SSSTT .* (dS.^3) .* (dT.^2);
        t_ssttt = (1/12) * R_SSTTT .* (dS.^2) .* (dT.^3);
        t_stttt = (1/24) * R_STTTT .* dS .* (dT.^4);
        
        calc_mm = @(rho_term) -(nansum(rho_term .* dz_perm, 3) / rho0) * 1000;
        
        TS_2_ST(:,:,idx) = calc_mm(t_st);
        TS_3_SST(:,:,idx) = calc_mm(t_sst); TS_3_STT(:,:,idx) = calc_mm(t_stt);
        TS_4_SSST(:,:,idx) = calc_mm(t_ssst); TS_4_SSTT(:,:,idx) = calc_mm(t_sstt); TS_4_STTT(:,:,idx) = calc_mm(t_sttt);
        TS_5_SSSST(:,:,idx) = calc_mm(t_sssst); TS_5_SSSTT(:,:,idx) = calc_mm(t_ssstt); TS_5_SSTTT(:,:,idx) = calc_mm(t_ssttt); TS_5_STTTT(:,:,idx) = calc_mm(t_stttt);
        
        if mod(idx, 24)==0, fprintf('   进度: %.1f%%\n', idx/TotalFiles*100); end
    end
end

%% 6. 统计趋势与振幅
fprintf('>> [4/7] 统计趋势与振幅...\n');

[Tr_ST, Am_ST] = calc_stats(TS_2_ST, time_vec, Nx, Ny);
[Tr_SST, Am_SST] = calc_stats(TS_3_SST, time_vec, Nx, Ny);
[Tr_STT, Am_STT] = calc_stats(TS_3_STT, time_vec, Nx, Ny);

[Tr_4_1, Am_4_1] = calc_stats(TS_4_SSST, time_vec, Nx, Ny);
[Tr_4_2, Am_4_2] = calc_stats(TS_4_SSTT, time_vec, Nx, Ny);
[Tr_4_3, Am_4_3] = calc_stats(TS_4_STTT, time_vec, Nx, Ny);

[Tr_5_1, Am_5_1] = calc_stats(TS_5_SSSST, time_vec, Nx, Ny);
[Tr_5_2, Am_5_2] = calc_stats(TS_5_SSSTT, time_vec, Nx, Ny);
[Tr_5_3, Am_5_3] = calc_stats(TS_5_SSTTT, time_vec, Nx, Ny);
[Tr_5_4, Am_5_4] = calc_stats(TS_5_STTTT, time_vec, Nx, Ny);

%% 7. 绘图准备
fprintf('>> [5/7] 准备绘图数据...\n');
Lon_Plot = Lon_Orig;
Lon_Plot(Lon_Plot > 180) = Lon_Plot(Lon_Plot > 180) - 360;
[Lon_Plot, sort_idx] = sort(Lon_Plot);
[Lon_Grid, Lat_Grid] = meshgrid(Lon_Plot, Lat);

% 经度重排的匿名函数（用于数据处理，不用于绘图命令串联）
reorder_func = @(x) x(:, sort_idx);

%% ================== 窗口 1: 低阶 (2-3阶) ==================
fprintf('>> [6/7] 绘制窗口 1 (2-3阶)...\n');
figure('Name', 'Low Order (2-3)', 'Position', [10, 50, 1400, 900], 'Color', 'w');
t1 = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

C_Div = m_colmap('diverging');
C_Amp = 'parula';

% 调用本地函数绘图（安全稳健）
% 1. ST
plot_single_panel(Tr_ST, '2阶趋势: \Delta S \Delta T', C_Div, [-0.5 0.5], 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_ST, '2阶振幅: \Delta S \Delta T', C_Amp, [0 5.0], 'mm', Lon_Grid, Lat_Grid, reorder_func);

% 2. SST
plot_single_panel(Tr_SST, '3阶趋势: 1/2 \Delta S^2 \Delta T', C_Div, [-0.2 0.2], 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_SST, '3阶振幅: 1/2 \Delta S^2 \Delta T', C_Amp, [0 0.5], 'mm', Lon_Grid, Lat_Grid, reorder_func);

% 3. STT
plot_single_panel(Tr_STT, '3阶趋势: 1/2 \Delta S \Delta T^2', C_Div, [-0.2 0.2], 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_STT, '3阶振幅: 1/2 \Delta S \Delta T^2', C_Amp, [0 0.5], 'mm', Lon_Grid, Lat_Grid, reorder_func);

sgtitle(t1, '低阶混合项 (2阶 & 3阶) 动力特征', 'FontSize', 16, 'FontWeight', 'bold');

%% ================== 窗口 2: 高阶 (4-5阶) ==================
fprintf('>> [7/7] 绘制窗口 2 (4-5阶)...\n');
figure('Name', 'High Order (4-5)', 'Position', [50, 50, 1800, 1000], 'Color', 'w');
t2 = tiledlayout(4, 4, 'TileSpacing', 'tight', 'Padding', 'compact');

% 智能获取色标范围
get_clim = @(x) [-1 1] * prctile(abs(x(~isnan(x))), 99);
get_alim = @(x) [0 1] * prctile(x(~isnan(x)), 99);

% --- Row 1: 4阶 趋势 ---
plot_single_panel(Tr_4_1, '4阶趋势: S^3 T', C_Div, get_clim(Tr_4_1), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Tr_4_2, '4阶趋势: S^2 T^2', C_Div, get_clim(Tr_4_2), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Tr_4_3, '4阶趋势: S T^3', C_Div, get_clim(Tr_4_3), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
nexttile; axis off; text(0,0.5,'4阶项\n(系数 1/6, 1/4, 1/6)','FontSize',14,'HorizontalAlignment','left');

% --- Row 2: 4阶 振幅 ---
plot_single_panel(Am_4_1, '4阶振幅: S^3 T', C_Amp, get_alim(Am_4_1), 'mm', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_4_2, '4阶振幅: S^2 T^2', C_Amp, get_alim(Am_4_2), 'mm', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_4_3, '4阶振幅: S T^3', C_Amp, get_alim(Am_4_3), 'mm', Lon_Grid, Lat_Grid, reorder_func);
nexttile; axis off; 

% --- Row 3: 5阶 趋势 ---
plot_single_panel(Tr_5_1, '5阶趋势: S^4 T', C_Div, get_clim(Tr_5_1), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Tr_5_2, '5阶趋势: S^3 T^2', C_Div, get_clim(Tr_5_2), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Tr_5_3, '5阶趋势: S^2 T^3', C_Div, get_clim(Tr_5_3), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Tr_5_4, '5阶趋势: S T^4', C_Div, get_clim(Tr_5_4), 'mm/yr', Lon_Grid, Lat_Grid, reorder_func);

% --- Row 4: 5阶 振幅 ---
plot_single_panel(Am_5_1, '5阶振幅: S^4 T', C_Amp, get_alim(Am_5_1), 'mm', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_5_2, '5阶振幅: S^3 T^2', C_Amp, get_alim(Am_5_2), 'mm', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_5_3, '5阶振幅: S^2 T^3', C_Amp, get_alim(Am_5_3), 'mm', Lon_Grid, Lat_Grid, reorder_func);
plot_single_panel(Am_5_4, '5阶振幅: S T^4', C_Amp, get_alim(Am_5_4), 'mm', Lon_Grid, Lat_Grid, reorder_func);

sgtitle(t2, '高阶混合项 (4阶 & 5阶) - 注意色标范围已自动缩小', 'FontSize', 16, 'FontWeight', 'bold');
fprintf('>> ✅ 所有绘图完成！请查看两个独立窗口。\n');

%% =======================================================
%% 本地函数定义 (全部放在文件最末尾)
%% =======================================================

% 1. 统计计算函数 (标准 for 循环)
function [Trend, Amp] = calc_stats(TS, time, Nx, Ny)
    Trend = nan(Ny, Nx); Amp = nan(Ny, Nx);
    for i = 1:Nx 
        for j = 1:Ny
            y = squeeze(TS(i,j,:));
            valid = ~isnan(y);
            if sum(valid) > length(time)*0.8
                p = polyfit(time(valid), y(valid), 1);
                Trend(j,i) = p(1);
                Amp(j,i) = std(y(valid) - polyval(p, time(valid)));
            end
        end
    end
end

% 2. 绘图单元函数 (解决 shading 报错的关键)
function plot_single_panel(data, title_str, cmap, clim_range, unit, Lon, Lat, reorder)
    nexttile;
    m_proj('robinson','lon',[-180 180],'lat',[-90 90]);
    m_pcolor(Lon, Lat, reorder(data));
    shading flat; % 这里是独立语句，完全合法
    m_coast('patch',[.9 .9 .9]);
    m_grid('linestyle',':','xticklabels',[],'yticklabels',[]);
    colormap(gca, cmap);
    caxis(clim_range);
    cb = colorbar;
    title(cb, unit);
    title(title_str, 'Interpreter', 'tex', 'FontSize', 12);
end