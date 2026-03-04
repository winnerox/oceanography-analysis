%% Analysis_Convergence_Loop_Avg.m
% =========================================================================
% 功能：平均态收敛性 Loop 分析 (Order 2 ~ 8)
% 比较：Exact vs Approx_k (k=2..8)
% 输出：
%   1. Trend Residual Maps (k=2..8)
%   2. Amp Residual Maps (k=2..8)
%   3. Convergence Boxplot
% 需要所有数据ready：T(1-8), S(1-8), Cross(Detailed)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
TermsDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\Convergence_Loop_Avg';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 加载数据
fprintf('>> [1/4] 加载各阶数据...\n');

% 2.1 Exact
load(fullfile(TermsDir, 'EN4_Formula11_Exact_Avg.mat'), 'TSLA_Exact_Avg', 'lon', 'lat', 'time_vec');
TSLA_Exact = TSLA_Exact_Avg;
Lon = lon; Lat = lat;
if max(Lon)>180, Lon(Lon>180)=Lon(Lon>180)-360; end
[Lon, sort_idx] = sort(Lon);
TSLA_Exact = TSLA_Exact(sort_idx, :, :);

% 2.2 T Terms (1-8)
load(fullfile(DataDir, 'EN4_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders');
TSLA_T = TSLA_AllOrders(sort_idx, :, :, 1:8);

% 2.3 S Terms (1-8) - Created by updated Calc_Cross_and_S_Terms_Avg.m
f_S = fullfile(TermsDir, 'EN4_Terms_S_Pure_Avg_1to8.mat');
if ~exist(f_S, 'file'), error('请先运行 Calc_Cross_and_S_Terms_Avg.m 生成 S项 1-8阶数据'); end
DS = load(f_S);

% 2.4 Cross Terms (Detailed) - Created by Calc_CrossDetail_Avg.m (MaxOrder=8 updated)
f_C = fullfile(TermsDir, 'EN4_CrossDetail_Avg.mat');
if ~exist(f_C, 'file'), error('请先运行 Calc_CrossDetail_Avg.m (MaxOrder=8) 生成混合项数据'); end
DC = load(f_C);

% Prepare Data Aggregation
[Nx, Ny, Nt] = size(TSLA_Exact);
[Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);

% Calc Function
compute_trend_amp = @(data, t) compute_stats_vectorized(data, t);
t_vec = time_vec(:);

%% 3. Loop Analysis (Order 2 to 8)
fprintf('>> [2/4] Loop Analysis (Orders 2-8)...\n');

Models = {'EN4', 'IAP', 'Ishii'};
% Models = {'EN4'}; % Debug
ModelColors = {'k', 'b', 'r'};

for m_idx = 1:length(Models)
    ModelName = Models{m_idx};
    fprintf('\n>> Processing Model: %s ...\n', ModelName);
    
    % --- Load Data ---
    % 1. Exact
    if strcmp(ModelName, 'EN4')
        fn = fullfile(TermsDir, 'EN4_Formula11_Exact_Avg.mat');
    else
        fn = fullfile(strrep(TermsDir, 'EN4', ModelName), sprintf('%s_Formula11_Exact_Avg.mat', ModelName));
    end
    if ~exist(fn, 'file'), fprintf('   [Skip] Missing Exact: %s\n', fn); continue; end
    D_Exact = load(fn);
    TSLA_Exact = D_Exact.TSLA_Exact_Avg; 
    Lon = D_Exact.lon; Lat = D_Exact.lat; Time = D_Exact.time_vec;
    
    % Fix Lon
    Lon = double(Lon); 
    Lon(Lon>180) = Lon(Lon>180)-360;
    [Lon, sort_idx] = sort(Lon);
    TSLA_Exact = TSLA_Exact(sort_idx, :, :);
    [Nx, Ny, Nt] = size(TSLA_Exact);
    [Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);
    t_vec = Time(:);

    % 2. T Terms (1-8)
    if strcmp(ModelName, 'EN4')
        % EN4 legacy file
        fn_T = fullfile(DataDir, 'EN4_TSLA_Terms1to8_Average.mat');
        if exist(fn_T, 'file')
            DT = load(fn_T, 'TSLA_AllOrders');
            TSLA_T_All = DT.TSLA_AllOrders(sort_idx, :, :, 1:8);
        else
             fprintf('   [Skip] Missing EN4 T Terms: %s\n', fn_T); continue;
        end
    else
        % IAP/Ishii New format
        fn_T = fullfile(strrep(TermsDir, 'EN4', ModelName), sprintf('%s_Terms_T_Pure_Avg_1to8.mat', ModelName));
        if ~exist(fn_T, 'file'), fprintf('   [Skip] Missing T Terms: %s\n', fn_T); continue; end
        DT = load(fn_T);
        % Convert cell to array for convenience or keep as cell? 
        % Code below expects TSLA_T(:,:,:,k)
        TSLA_T_All = zeros(Nx, Ny, Nt, 8);
        for k=1:8
            vn = sprintf('TTLA_%d', k);
            if isfield(DT, vn)
                TSLA_T_All(:,:,:,k) = DT.(vn)(sort_idx, :, :);
            end
        end
    end

    % 3. S Terms (1-8)
    if strcmp(ModelName, 'EN4')
        fn_S = fullfile(TermsDir, 'EN4_Terms_S_Pure_Avg_1to8.mat');
    else
        fn_S = fullfile(strrep(TermsDir, 'EN4', ModelName), sprintf('%s_Terms_S_Pure_Avg_1to8.mat', ModelName));
    end
    if ~exist(fn_S, 'file'), fprintf('   [Skip] Missing S Terms: %s\n', fn_S); continue; end
    DS = load(fn_S);

    % 4. Cross Terms
    if strcmp(ModelName, 'EN4')
        fn_C = fullfile(TermsDir, 'EN4_CrossDetail_Avg.mat');
    else
        fn_C = fullfile(strrep(TermsDir, 'EN4', ModelName), sprintf('%s_CrossDetail_Avg.mat', ModelName));
    end
    if ~exist(fn_C, 'file'), fprintf('   [Skip] Missing Cross Terms: %s\n', fn_C); continue; end
    DC = load(fn_C);

    % --- Analyze Orders ---
    Trend_Res_All = cell(8,1);
    Amp_Res_All   = cell(8,1);
    
    Sum_T = zeros(Nx, Ny, Nt);
    Sum_S = zeros(Nx, Ny, Nt);
    Sum_Cross = zeros(Nx, Ny, Nt);

    for k = 1:8
        % Add T_k
        Sum_T = Sum_T + squeeze(TSLA_T_All(:,:,:,k));
        
        % Add S_k
        vn_s = sprintf('SSLA_%d', k);
        if isfield(DS, vn_s)
            Sum_S = Sum_S + DS.(vn_s)(sort_idx, :, :);
        end
        
        % Add Cross Terms (Order k: n+m=k)
        for n = 1:k-1
            m = k - n;
            vn_c = sprintf('Cross_T%dS%d', n, m);
            if isfield(DC, vn_c)
                Sum_Cross = Sum_Cross + DC.(vn_c)(sort_idx, :, :);
            end
        end
        
        % Approx & Residual
        Approx_k = Sum_T + Sum_S + Sum_Cross;
        Residual_k = TSLA_Exact - Approx_k;
        
        % Stats
        [tr, am] = compute_trend_amp(Residual_k, t_vec);
        Trend_Res_All{k} = tr;
        Amp_Res_All{k}   = am;
    end
    
    % --- Plotting per Model ---
    PlotOrders = 2:8; N_Plot = length(PlotOrders);
    
    % Trend Map
    figure('Position', [50, 50, 1600, 800], 'Color', 'w');
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:N_Plot
        ord = PlotOrders(i);
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        m_pcolor(Lon_Grid, Lat_Grid, Trend_Res_All{ord}); shading interp;
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4]);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.7 .7 .7], 'fontsize', 7, 'xticklabels',[],'yticklabels',[]);
        colormap(gca, MyDivergentCmap);
        cl = prctile(abs(Trend_Res_All{ord}(:)), 98); if cl<0.01, cl=0.01; end
        caxis([-cl, cl]);
        title(sprintf('Res Trend (Order %d)', ord), 'FontSize', 11);
        cb=colorbar; cb.Label.String='mm/yr';
    end
    sgtitle(sprintf('%s Avg State Residual Trends', ModelName), 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(OutputDir, sprintf('Res_Trend_Maps_%s.png', ModelName)));
    
    % Amp Map
    figure('Position', [100, 100, 1600, 800], 'Color', 'w');
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:N_Plot
        ord = PlotOrders(i);
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        m_pcolor(Lon_Grid, Lat_Grid, Amp_Res_All{ord}); shading interp;
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4]);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.7 .7 .7], 'fontsize', 7, 'xticklabels',[],'yticklabels',[]);
        colormap(gca, MySeqCmap);
        cl = prctile(Amp_Res_All{ord}(:), 98); if cl<0.01, cl=0.01; end
        caxis([0, cl]);
        title(sprintf('Res Amp (Order %d)', ord), 'FontSize', 11);
        cb=colorbar; cb.Label.String='mm';
    end
    sgtitle(sprintf('%s Avg State Residual Amplitudes', ModelName), 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(OutputDir, sprintf('Res_Amp_Maps_%s.png', ModelName)));
    
    % Boxplot (Individual Model)
    % ... (Optional, or aggregate later)
end

%% 5. Boxplot Summary
% Using draw_jitter_unit
fprintf('>> [4/4] 绘制 Boxplots...\n');
figure('Position', [200, 200, 1000, 500], 'Color', 'w');
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Trend
nexttile; hold on;
d1 = Trend_Res_All{2}(:); d8 = Trend_Res_All{8}(:);
ymax = max(abs([d1; d8; prctile(d1,99)*1.2]));
if ymax==0, ymax=1; end

for i = 1:N_Plot
    ord = PlotOrders(i);
    data = Trend_Res_All{ord}(:);
    % Color gradient from Red to Blue? Or just Blue
    col = [0.2 0.4 0.8];
    draw_jitter_unit(gca, ord, data, 0.5, 0.3, 10, col, [0.95 0.95 1], 0.3);
end
yline(0, 'k--');
xlim([1.5, 8.5]); xticks(2:8); xlabel('Order'); ylabel('Residual Trend (mm/yr)');
title('Convergence of Residual Trend');
ylim([-ymax, ymax]);

% Amp
nexttile; hold on;
ymax = max([Amp_Res_All{2}(:); prctile(Amp_Res_All{2}(:),99)*1.2]);
for i = 1:N_Plot
    ord = PlotOrders(i);
    data = Amp_Res_All{ord}(:);
    draw_jitter_unit(gca, ord, data, 0.5, 0.3, 10, [0.8 0.3 0.2], [1 0.95 0.95], 0.3);
end
xlim([1.5, 8.5]); xticks(2:8); xlabel('Order'); ylabel('Residual Amp (mm)');
title('Convergence of Residual Amplitude');
ylim([0, ymax]);

saveas(gcf, fullfile(OutputDir, 'Res_Boxplot_Loop.png'));

fprintf('>> Done.\n');

%% Helper
function [trend, amp] = compute_stats_vectorized(data_3d, t_vec)
    [Nx, Ny, Nt] = size(data_3d);
    t_vec = t_vec(:); 
    trend = nan(Nx, Ny); amp = nan(Nx, Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            y = squeeze(data_3d(i, j, :));
            mask = ~isnan(y);
            if sum(mask) < Nt * 0.7, continue; end
            p = polyfit(t_vec(mask), y(mask), 1);
            trend(i, j) = p(1);
            resid = y(mask) - polyval(p, t_vec(mask));
            amp(i, j) = std(resid);
        end
    end
    trend = trend'; amp = amp';
end
