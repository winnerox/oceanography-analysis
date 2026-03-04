%% Analysis_Convergence_Loop_Std.m
% =========================================================================
% 功能：标准态收敛性 Loop 分析 (Order 2 ~ 8)
% 比较：Exact vs Approx_k (k=2..8)
% =========================================================================
clear; clc; close all;

%% 1. 配置
DataDir = 'D:\work\MAT_Data';
TermsDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\Figures\Convergence_Loop_Std';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Models = {'EN4', 'IAP', 'Ishii'};
% Models = {'EN4'}; 
ModelColors = {'k', 'b', 'r'};

%% 2. Loop Models
for m_idx = 1:length(Models)
    ModelName = Models{m_idx};
    fprintf('\n>> Processing Model: %s ...\n', ModelName);
    
    % --- Determine Paths ---
    if strcmp(ModelName, 'EN4')
        ModelTermsDir = TermsDir; % D:\work\EN4_TSLA_Terms
        fn_Exact = fullfile(ModelTermsDir, 'EN4_Formula11_Exact_Std.mat');
        fn_T     = fullfile(DataDir, 'EN4_TSLA_Terms_1to8_StdRef.mat');
        fn_Bundled = fullfile(ModelTermsDir, 'EN4_Std_Terms_HighOrder_1to8.mat'); % S + Cross
    else
        ModelTermsDir = strrep(TermsDir, 'EN4', ModelName); % D:\work\IAP_TSLA_Terms
        fn_Exact = fullfile(ModelTermsDir, sprintf('%s_Formula11_Exact_Std.mat', ModelName));
        fn_T     = fullfile(ModelTermsDir, sprintf('%s_Terms_T_Pure_Std_1to8.mat', ModelName));
        fn_S     = fullfile(ModelTermsDir, sprintf('%s_Terms_S_Pure_Std_1to8.mat', ModelName));
        fn_C     = fullfile(ModelTermsDir, sprintf('%s_CrossDetail_Std.mat', ModelName));
    end
    
    if ~exist(fn_Exact, 'file'), fprintf('   [Skip] Missing Exact: %s\n', fn_Exact); continue; end
    
    % --- Load Data ---
    % 1. Exact
    D_Exact = load(fn_Exact);
    if isfield(D_Exact, 'TSLA_Exact_Std')
        TSLA_Exact = D_Exact.TSLA_Exact_Std;
    else
        % Try dynamic field name? Or assume fixed.
        TSLA_Exact = D_Exact.TSLA_Exact_Std; 
    end
    Lon = D_Exact.lon; Lat = D_Exact.lat; Time = D_Exact.time_vec;
    
    % Fix Lon Logic for Analysis (Uniform -180~180)
    Lon = double(Lon); 
    Lon(Lon>180) = Lon(Lon>180)-360;
    [Lon, sort_idx] = sort(Lon);
    TSLA_Exact = TSLA_Exact(sort_idx, :, :);
    
    [Nx, Ny, Nt] = size(TSLA_Exact);
    [Lon_Grid, Lat_Grid] = meshgrid(Lon, Lat);
    t_vec = Time(:);
    compute_trend_amp = @(data, t) compute_stats_vectorized(data, t);

    % 2. T Terms
    if ~exist(fn_T, 'file'), fprintf('   [Skip] Missing T Terms: %s\n', fn_T); continue; end
    DT = load(fn_T);
    TSLA_T_All = zeros(Nx, Ny, Nt, 8);
    
    if strcmp(ModelName, 'EN4')
         % EN4 has 4D array
         TSLA_T_All = DT.TSLA_AllOrders(sort_idx, :, :, 1:8);
    else
         % IAP/Ishii has TTLA_k
         for k=1:8
            vn = sprintf('TTLA_%d', k);
            if isfield(DT, vn)
                TSLA_T_All(:,:,:,k) = DT.(vn)(sort_idx, :, :);
            end
         end
    end
    
    % 3. S and Cross Terms
    ds_s = []; ds_c = [];
    if strcmp(ModelName, 'EN4')
        if ~exist(fn_Bundled, 'file'), fprintf('   [Skip] Missing Bundled SC Terms: %s\n', fn_Bundled); continue; end
        ds_bundled = load(fn_Bundled);
        ds_s = ds_bundled; 
        ds_c = ds_bundled; % Both in same struct
    else
        if ~exist(fn_S, 'file'), fprintf('   [Skip] Missing S Terms: %s\n', fn_S); continue; end
        if ~exist(fn_C, 'file'), fprintf('   [Skip] Missing Cross Terms: %s\n', fn_C); continue; end
        ds_s = load(fn_S);
        ds_c = load(fn_C);
    end
    
    % --- Loop Analysis ---
    Trend_Res_All = cell(8,1);
    Amp_Res_All   = cell(8,1);
    Sum_Total = zeros(Nx, Ny, Nt); % Cumulative sum of Approx
    
    fprintf('   Calcing Orders 2..8: ');
    for k = 1:8
        term_k = squeeze(TSLA_T_All(:,:,:,k));
        
        % Add S Term
        vn_s = sprintf('SSLA_%d', k);
        if isfield(ds_s, vn_s)
            term_k = term_k + ds_s.(vn_s)(sort_idx, :, :);
        end
        
        % Add Cross Term (Bundled logic for EN4 uses Cross_k, Split uses Cross_TnSm)
        if strcmp(ModelName, 'EN4')
             if k>=2
                vn_c = sprintf('Cross_%d', k);
                if isfield(ds_c, vn_c)
                    term_k = term_k + ds_c.(vn_c)(sort_idx, :, :);
                end
             end
        else
            % IAP/Ishii (Sum all combinations for order k)
            if k>=2
                for n = 1:k-1
                    m = k - n;
                    vn_c = sprintf('Cross_T%dS%d', n, m);
                    if isfield(ds_c, vn_c)
                        term_k = term_k + ds_c.(vn_c)(sort_idx, :, :);
                    end
                end
            end
        end
        
        Sum_Total = Sum_Total + term_k;
        resid = TSLA_Exact - Sum_Total;
        
        [tr, am] = compute_trend_amp(resid, t_vec);
        Trend_Res_All{k} = tr;
        Amp_Res_All{k}   = am;
    end
    fprintf('Done.\n');
    
    % --- Plotting per Model ---
    PlotOrders = 2:8; N_Plot = length(PlotOrders);
    
    % Helper colormaps
    MyBlue=[0.05 0.25 0.6]; MyWhite=[1 1 1]; MyRed=[0.7 0.05 0.05];
    xx=linspace(0,1,256);
    MyDivergentCmap = [interp1([0 0.5 1], [MyBlue(1); MyWhite(1); MyRed(1)], xx)', ...
                       interp1([0 0.5 1], [MyBlue(2); MyWhite(2); MyRed(2)], xx)', ...
                       interp1([0 0.5 1], [MyBlue(3); MyWhite(3); MyRed(3)], xx)'];
    MyRedSeq = [1 1 1; 1 0.6 0.4; 0.4 0 0];
    MySeqCmap = [interp1([0 0.5 1], MyRedSeq(:,1), xx)', interp1([0 0.5 1], MyRedSeq(:,2), xx)', interp1([0 0.5 1], MyRedSeq(:,3), xx)'];

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
    sgtitle(sprintf('%s Std State Residual Trends', ModelName), 'FontSize', 14, 'FontWeight', 'bold');
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
    sgtitle(sprintf('%s Std State Residual Amplitudes', ModelName), 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(OutputDir, sprintf('Res_Amp_Maps_%s.png', ModelName)));
    
    drawnow;
end

fprintf('>> Done.\n');

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
