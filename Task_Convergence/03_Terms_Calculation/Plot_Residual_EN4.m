%% Plot_Residual_EN4.m
% 混合项前 N 阶截断近似与精确海平面的残差绘图 (按照参考图样式)
% 模型: EN4
% 输出: Map_Amp, Map_Trend, Boxplot_Amp, Boxplot_Trend 的 2~8 阶残差
clear; clc; close all;
ws = warning('off', 'all');
addpath('D:\work');
addpath('D:\work\m_map');

%% 1. 配置参数
Model = 'EN4';
DataDir = 'D:\work\EN4_TSLA_Terms';
OutputDir = 'D:\work\MAT_Data\Figures\Residuals';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

States = {'Avg', 'Std'};
ExactFileSuffix = {'Formula11_Exact_Avg.mat', 'Formula11_Exact_Std.mat'};
TermsFileSuffix = {'TSLA_Terms_1to8_Average.mat', 'TSLA_Terms_1to8_StdRef.mat'};
ExactVarName = {'TSLA_Exact_Avg', 'TSLA_Exact_Std'};

%% 配色定义
MyRed = [1 1 1; 1 0.9 0.8; 1 0.6 0.4; 0.9 0.1 0.1; 0.4 0 0];
xx = [0 0.1 0.4 0.8 1]; yy = linspace(0,1,256);
MySeqCmap = [interp1(xx,MyRed(:,1),yy)', interp1(xx,MyRed(:,2),yy)', interp1(xx,MyRed(:,3),yy)'];

C_Blue = [0.05 0.25 0.6];   C_LBlue = [0.6 0.8 1]; C_White = [1 1 1];
C_LRed = [1 0.7 0.6];       C_Red2 = [0.7 0.05 0.05];
Locs = [0, 0.25, 0.5, 0.75, 1]; 
CPts = [C_Blue; C_LBlue; C_White; C_LRed; C_Red2];
xx2 = linspace(0, 1, 256);
MyDivCmap = [interp1(Locs, CPts(:,1), xx2, 'pchip')', ...
             interp1(Locs, CPts(:,2), xx2, 'pchip')', ...
             interp1(Locs, CPts(:,3), xx2, 'pchip')'];

BoxColors = {[0.85 0.325 0.098], [0 0.447 0.741], [0.466 0.674 0.188]};
BoxBgColor = [0.95 0.95 1];

for s_idx = 1:2
    state = States{s_idx};
    fprintf('\n============= 处理 %s (%s) =============\n', Model, state);
    
    ExactFile = fullfile(DataDir, sprintf('%s_%s', Model, ExactFileSuffix{s_idx}));
    TermsFile = fullfile(DataDir, sprintf('%s_%s', Model, TermsFileSuffix{s_idx}));
    if ~exist(ExactFile, 'file') || ~exist(TermsFile, 'file')
        fprintf('缺少 %s 数据文件, 跳过.\n', state); continue;
    end
    
    %% 加载与预处理
    tmp_e = load(ExactFile, ExactVarName{s_idx});
    ExactData = single(tmp_e.(ExactVarName{s_idx}));
    
    % 自动识别变量名，因为 Avg 和 Std 保存的格式可能有大小写区别
    tmp_all = load(TermsFile);
    
    if isfield(tmp_all, 'Lon'), lon_v = tmp_all.Lon; else, lon_v = tmp_all.lon; end
    if isfield(tmp_all, 'Lat'), lat_v = tmp_all.Lat; else, lat_v = tmp_all.lat; end
    if isfield(tmp_all, 'Time_Axis'), time_v = tmp_all.Time_Axis; 
    elseif isfield(tmp_all, 'time_axis'), time_v = tmp_all.time_axis;
    elseif isfield(tmp_all, 'time_vec'), time_v = tmp_all.time_vec; 
    else, time_v = []; end
    
    if isfield(tmp_all, 'TSLA_AllOrders')
        terms_v = tmp_all.TSLA_AllOrders;
    elseif isfield(tmp_all, 'TSLA_Terms_AllOrders')
        terms_v = tmp_all.TSLA_Terms_AllOrders;
    end
    
    lon = wrapTo180(lon_v); lat = lat_v; t_vec = time_v(:);
    [lon, sort_idx] = sort(lon);
    [Lon_Grid, Lat_Grid] = meshgrid(lon, lat);
    
    TermsData = single(terms_v(sort_idx, :, :, :));
    ExactData = ExactData(sort_idx, :, :);
    
    [Nx, Ny, Nt, ~] = size(TermsData);
    
    %% 时间序列对齐检查 (自动纠正可能存在的时间偏移)
    % 采用对全局空间求均值获得1D时间序列，避免空间索引不一致造成的点错位
    fprintf('>> 检查时间序列对齐...\n');
    ts_e_global = squeeze(nanmean(ExactData, [1, 2]));
    ts_t_global = squeeze(nanmean(TermsData(:,:,:,1), [1, 2]));
    
    valid_mask = ~isnan(ts_e_global) & ~isnan(ts_t_global);
    if sum(valid_mask) > 50
        ts_e_clean = ts_e_global(valid_mask); ts_e_clean = ts_e_clean - nanmean(ts_e_clean); 
        ts_t_clean = ts_t_global(valid_mask); ts_t_clean = ts_t_clean - nanmean(ts_t_clean);
        [xc, lags] = xcorr(ts_e_clean, ts_t_clean, 10);
        [~, max_idx] = max(xc);
        optimal_lag = lags(max_idx);
        
        if optimal_lag ~= 0
            fprintf('⚠️ 检测到 ExactData 与 TermsData 存在 %d 个月的时间偏移！执行自动对齐...\n', optimal_lag);
            if optimal_lag > 0
                ExactData = ExactData(:,:, optimal_lag+1:end);
                TermsData = TermsData(:,:, 1:end-optimal_lag, :);
                t_vec = t_vec(optimal_lag+1:end);
            else
                ExactData = ExactData(:,:, 1:end+optimal_lag);
                TermsData = TermsData(:,:, -optimal_lag+1:end, :);
                t_vec = t_vec(1:end+optimal_lag);
            end
            [Nx, Ny, Nt, ~] = size(TermsData);
        else
            fprintf('✅ 时间序列已原生对齐.\n');
        end
    else
        fprintf('⚠️ 可用于时间序列对齐的有效点太少.\n');
    end
    
    %% 空间网格完美对齐 (解决 Nx, Ny, 倒换或部分边界被 NaN 掩盖导致相减完全错位的问题)
    % ExactData 此时尺寸为 [N_lon_e, N_lat_e, Nt]，TermsData 为 [N_lon_t, N_lat_t, Nt] (前面可能已经被 sort_idx 倒换)
    % 为了确保万无一失，直接用 interp2 把 TermsData 插值到 ExactData 的网格上
    fprintf('>> 检查空间网格对齐...\n');
    [Ny_e, Nx_e, Nt_e] = size(ExactData);
    
    if ~isequal(size(ExactData), size(TermsData(:,:,:,1))) || any(lon_v ~= lon_v) % 强制检查
        fprintf('⚠️ 检测到空间维度排序不一致，执行严格的网格重采样对齐...\n');
        
        Aligned_TermsData = nan(Ny_e, Nx_e, Nt_e, 8, 'single');
        [X_e, Y_e] = ndgrid(lon, lat); % ExactGrid (sort过的lon)
        [X_t, Y_t] = ndgrid(lon_v, lat_v); % Terms 原版网格
        
        % 如果 TermsData 前面已经被应用了 sort_idx，现在用原版再倒一次会出错，所以用原版数据
        tmp_terms_raw = single(terms_v); 
        
        for n_order = 1:8
            for t_idx = 1:Nt_e
                % 提取当前切片并在原始网格上插值到 Exact 目标网格
                slice = squeeze(tmp_terms_raw(:,:,t_idx,n_order));
                F = griddedInterpolant(X_t, Y_t, slice, 'linear', 'none');
                Aligned_TermsData(:,:,t_idx,n_order) = F(X_e, Y_e);
            end
        end
        TermsData = Aligned_TermsData;
    end
    [Nx, Ny, Nt, ~] = size(TermsData);
    
    %% 单位检测 (基于原始信号，而不是残差，避免将收敛后的极小值误扩 1000 倍)
    check_a = nanmean(std(ExactData, 0, 3, 'omitnan'), 'all');
    need_scale = (check_a < 0.5);
    if need_scale
        ExactData = ExactData * 1000;
        TermsData = TermsData * 1000;
    end
    
    %% 计算 1~8 阶累加近似与残差
    Residual_Amp = nan(Ny, Nx, 8);
    Residual_Trend = nan(Ny, Nx, 8);
    Residual_Sig = nan(Ny, Nx, 8);
    
    valid_mask = ~isnan(ExactData) & ~all(isnan(TermsData), 4) & (ExactData ~= 0);
    Sum_Taylor = zeros(Nx, Ny, Nt, 'single');
    
    fprintf('>> 计算残差时空序列 (1-8阶)...\n');
    Data_Amp_Box = cell(8, 1);
    Data_Trend_Box = cell(8, 1);
    
    for n = 1:8
        term_n = single(TermsData(:,:,:,n)); term_n(isnan(term_n)) = 0;
        Sum_Taylor = Sum_Taylor + term_n;
        
        Res_3D = ExactData - Sum_Taylor;
        Res_3D(~valid_mask) = NaN;
        
        temp_amp = nan(Nx, Ny);
        temp_trend = nan(Nx, Ny);
        temp_sig = nan(Nx, Ny);
        
        box_amp_list = [];
        box_trend_list = [];
        
        for i = 1:Nx
            for j = 1:Ny
                ts = squeeze(Res_3D(i,j,:));
                idx_v = ~isnan(ts);
                if sum(idx_v) < Nt * 0.5, continue; end
                
                [p_coeff, ~] = polyfit(t_vec(idx_v), ts(idx_v), 1);
                y_fit = polyval(p_coeff, t_vec(idx_v));
                detrended_ts = ts(idx_v) - y_fit;
                
                ts_amp = std(detrended_ts);
                ts_trend = p_coeff(1);
                
                temp_amp(i,j) = ts_amp;
                temp_trend(i,j) = ts_trend;
                
                box_amp_list(end+1) = ts_amp; %#ok<AGROW>
                box_trend_list(end+1) = ts_trend; %#ok<AGROW>
                
                % t-test
                resid = ts(idx_v) - y_fit; n_val = sum(idx_v);
                X = [ones(n_val,1) t_vec(idx_v)];
                XTX_inv = inv(X' * X);
                se_slope = sqrt(var(resid) * XTX_inv(2,2));
                t_stat = p_coeff(1) / se_slope;
                p_val = 2 * (1 - tcdf(abs(t_stat), n_val - 2));
                if p_val < 0.05, temp_sig(i,j) = 1; end
            end
        end
        Residual_Amp(:,:,n) = temp_amp';
        Residual_Trend(:,:,n) = temp_trend';
        Residual_Sig(:,:,n) = temp_sig';
        
        Data_Amp_Box{n} = box_amp_list;
        Data_Trend_Box{n} = box_trend_list;
        fprintf('  完成阶数: %d\n', n);
    end
    
    % (原本依赖残差幅度缩放 1000 倍的错误逻辑已被移除，改在了数据相减前对原始信号进行统一缩放)
    
    % 只画 2~8 阶 (对应剩余误差)，共7张子图，我们凑个 4x2 布局，留空一格或者第一格放 exact
    PlotOrders = 1:8; % 按照源模板习惯画8张
    
    %% 出图 1: Map_Amp
    f1 = figure('Visible', 'on', 'Position', [100, 50, 1200, 1400], 'Color', 'w');
    t1 = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for n = PlotOrders
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        data = Residual_Amp(:,:,n);
        m_pcolor(Lon_Grid, Lat_Grid, data); shading flat;
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
               'linewidth', 0.5, 'fontsize', 8);
        colormap(gca, MySeqCmap);
        
        raw_max = prctile(data(:), 95.0); if isnan(raw_max)||raw_max<1e-12, raw_max=1e-12; end
        expo = floor(log10(raw_max)); frac = raw_max / 10^expo;
        if frac <= 1
            nice_base = 1;
        elseif frac <= 1.25
            nice_base = 1.25;
        elseif frac <= 1.5
            nice_base = 1.5;
        elseif frac <= 2
            nice_base = 2;
        elseif frac <= 2.5
            nice_base = 2.5;
        elseif frac <= 3
            nice_base = 3;
        elseif frac <= 4
            nice_base = 4;
        elseif frac <= 5
            nice_base = 5;
        elseif frac <= 6
            nice_base = 6;
        elseif frac <= 8
            nice_base = 8;
        else
            nice_base = 10;
        end
        flim = nice_base * 10^expo;
        if flim == 0, flim = 1; end
        caxis([0, flim]);
        
        title(sprintf('Residual Amp - Sum(1..%d)', n), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        cb = colorbar; cb.FontSize = 8; 
        title(cb, 'mm', 'FontSize', 8, 'FontName', 'Times New Roman');
        cb.Ticks = linspace(0, flim, 5);
        try cb.Ruler.Exponent = 0; catch; end
    end
    sgtitle(sprintf('\\fontname{Microsoft YaHei}%s (%s) 累加项与精确解振幅空间残差图', Model, state), 'FontSize', 16, 'FontWeight', 'bold');
    set(findall(f1, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(f1, 'PaperPositionMode', 'auto'); 
    exportgraphics(f1, fullfile(OutputDir, sprintf('%s_%s_Map_Amp_Res.png', Model, state)));
    % close(f1); % 注释掉以免立刻关闭图窗
    
    %% 出图 2: Map_Trend
    f2 = figure('Visible', 'on', 'Position', [100, 50, 1200, 1400], 'Color', 'w');
    t2 = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for n = PlotOrders
        nexttile;
        m_proj('robinson', 'lon', [-180 180], 'lat', [-90 90]);
        data = Residual_Trend(:,:,n); sig = Residual_Sig(:,:,n);
        m_pcolor(Lon_Grid, Lat_Grid, data); shading flat; hold on;
        
        skip = 3; [row, ~] = find(sig(1:skip:end, 1:skip:end) == 1);
        if ~isempty(row)
            slon = Lon_Grid(1:skip:end, 1:skip:end); slat = Lat_Grid(1:skip:end, 1:skip:end);
            m_plot(slon(sig(1:skip:end,1:skip:end)==1), slat(sig(1:skip:end,1:skip:end)==1), ...
                'o', 'markerfacecolor', [.2 .2 .2], 'markeredgecolor', 'none', 'markersize', 0.8);
        end
        m_coast('patch', [.92 .92 .92], 'edgecolor', [.4 .4 .4], 'linewidth', 0.5);
        m_grid('tickdir', 'out', 'linest', ':', 'gridcolor', [.65 .65 .65], ...
               'linewidth', 0.5, 'fontsize', 8);
        colormap(gca, MyDivCmap);
        
        raw_max = prctile(abs(data(:)), 99.8); if isnan(raw_max)||raw_max<1e-12, raw_max=1e-12; end
        expo = floor(log10(raw_max)); frac = raw_max / 10^expo;
        if frac<=2, bas=2; elseif frac<=5, bas=5; else, bas=10; end
        flim = bas * 10^expo;
        if flim == 0, flim = 1; end
        caxis([-flim, flim]);
        
        title(sprintf('Residual Trend - Sum(1..%d)', n), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        cb = colorbar; cb.FontSize = 8; 
        title(cb, 'mm/yr', 'FontSize', 8, 'FontName', 'Times New Roman');
        cb.Ticks = linspace(-flim, flim, 5);
        try cb.Ruler.Exponent = 0; catch; end
    end
    sgtitle(sprintf('\\fontname{Microsoft YaHei}%s (%s) 累加项与精确解趋势空间残差图 (打点 p<0.05)', Model, state), 'FontSize', 16, 'FontWeight', 'bold');
    set(findall(f2, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(f2, 'PaperPositionMode', 'auto'); 
    exportgraphics(f2, fullfile(OutputDir, sprintf('%s_%s_Map_Trend_Res.png', Model, state)));
    % close(f2); % 注释掉以免立刻关闭图窗
    
    %% 箱线图辅助函数
    draw_jit = @(ax, center, data, clr, cbox) ...
        draw_jitter(ax, center, data, clr, cbox); 

    %% 出图 3: Boxplot_Amp
    f3 = figure('Visible', 'on', 'Position', [50, 500, 1200, 550], 'Color', 'w');
    ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');
    for i = 1:8
        draw_jit(ax, i, Data_Amp_Box{i}, [0.2 0.4 0.7], BoxBgColor);
    end
    yline(ax, 0, '-', 'Color', 'k', 'LineWidth', 1.5);
    set(ax, 'Box', 'off', 'XTick', 1:8, 'XTickLabel', sprintfc('R_{%d}', 1:8));
    title(sprintf('%s (%s) 振幅残差序列分布 (单轴)', Model, state), 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Residual Amplitude (mm)');
    
    % 智能缩放上限
    all_d = horzcat(Data_Amp_Box{:});
    if ~isempty(all_d), ylim([0, max(prctile(all_d, 99.5), 1)*1.2]); end
    
    exportgraphics(f3, fullfile(OutputDir, sprintf('%s_%s_Boxplot_Amp_Res.png', Model, state)));
    % close(f3); % 注释掉以免立刻关闭图窗
    
    %% 出图 4: Boxplot_Trend
    f4 = figure('Visible', 'on', 'Position', [50, 50, 1200, 550], 'Color', 'w');
    ax = axes('Position', [0.07, 0.14, 0.90, 0.80]); hold(ax, 'on');
    for i = 1:8
        draw_jit(ax, i, Data_Trend_Box{i}, [0.2 0.4 0.7], BoxBgColor);
    end
    yline(ax, 0, '--', 'Color', [.5 .5 .5], 'LineWidth', 1.5);
    set(ax, 'Box', 'off', 'XTick', 1:8, 'XTickLabel', sprintfc('R_{%d}', 1:8));
    title(sprintf('%s (%s) 趋势残差序列分布 (单轴)', Model, state), 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Residual Trend (mm/yr)');
    
    all_d2 = horzcat(Data_Trend_Box{:});
    if ~isempty(all_d2), maxv = max(abs(prctile(all_d2, [0.5 99.5]))); ylim([-maxv, maxv]*1.2); end
    
    exportgraphics(f4, fullfile(OutputDir, sprintf('%s_%s_Boxplot_Trend_Res.png', Model, state)));
    % close(f4); % 注释掉以免立刻关闭图窗
    
end
fprintf('====== EN4 全部完成 ======\n');

function draw_jitter(ax, center, data, color_dots, color_box)
    data = data(~isnan(data)); if isempty(data), return; end
    q1 = prctile(data, 25); q3 = prctile(data, 75); med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w = prctile(data, 97.5);
    out = data(data < low_w | data > up_w);
    if length(out) > 2000, out = out(randperm(length(out), 2000)); end
    if ~isempty(out)
        xj = center + (rand(size(out)) - 0.5) * 0.4;
        scatter(ax, xj, out, 12, color_dots, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2);
    end
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    patch(ax, center+[-0.275 0.275 0.275 -0.275], [q1 q1 q3 q3], color_box, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot(ax, [center-0.275, center+0.275], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end
