%% Plot_Convergence_Residual_Map.m
% 绘制三家数据（EN4 / IAP / Ishii）在平均态/标准态下的 2~8阶 残差 RMSE 空间分布图
% 残差 = Exact_TSLA - sum(TSLA_Terms_AllOrders, 4)
clear; clc;
ws = warning('off', 'all');
addpath('D:\work');
addpath('D:\work\m_map'); % 假设存在 m_map

%% 配置
OutputFigureDir = 'D:\work\MAT_Data\Figures';
if ~exist(OutputFigureDir, 'dir'), mkdir(OutputFigureDir); end

Models = {'EN4', 'IAP', 'Ishii'};
States = {'Avg', 'Std'};
ExactFileSuffix = {'Formula11_Exact_Avg.mat', 'Formula11_Exact_Std.mat'};
TermsFileSuffix = {'TSLA_Terms_1to8_Average.mat', 'TSLA_Terms_1to8_StdRef.mat'};

Orders2Plot = 2:8; % 我们需要绘制的组合阶数截断点
N_Orders = length(Orders2Plot);

% 色标设定 (RMSE，从0到某个极大值)
caxis_range = [0, 5]; % mm (可以根据实际数据进行调整，保证不同图对比一致)
cmap = turbo(64);

%% 画图循环：分为 Avg 和 Std 两张大图
for s_idx = 1:length(States)
    state = States{s_idx};
    fprintf('正在处理 %s 状态残差图...\n', state);
    
    % 创建不可见的大图版 20寸宽 x 30寸高
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0 0 15 20], 'Color', 'w');
    
    for m_idx = 1:length(Models)
        model = Models{m_idx};
        
        DirName = sprintf('D:\\work\\%s_TSLA_Terms', model);
        ExactFile = fullfile(DirName, sprintf('%s_%s', model, ExactFileSuffix{s_idx}));
        TermsFile = fullfile(DirName, sprintf('%s_%s', model, TermsFileSuffix{s_idx}));
        
        if ~exist(ExactFile, 'file') || ~exist(TermsFile, 'file')
            fprintf('  警告: 缺少 %s %s 数据文件，跳过该列。\n', model, state);
            continue;
        end
        
        % 1. 读取 Exact (Lon x Lat x Time)
        fprintf('  读取 %s ...\n', model);
        tmp_e = load(ExactFile, 'TSLA_Exact_Avg');  % IAP和Ishii, EN4 变量名可能叫 TSLA_Exact_Avg
        if isfield(tmp_e, 'TSLA_Exact_Avg')
            Exact_TSLA = single(tmp_e.TSLA_Exact_Avg);
        else
            % 尝试读取可能的 Std 变量名
            tmp_s = load(ExactFile, 'TSLA_Exact_Std');
            if isfield(tmp_s, 'TSLA_Exact_Std')
                Exact_TSLA = single(tmp_s.TSLA_Exact_Std);
            else
                % 直接读第一个变量
                fns = fieldnames(tmp_e);
                for i = 1:length(fns)
                    if contains(fns{i}, 'TSLA') && isnumeric(tmp_e.(fns{i}))
                        Exact_TSLA = single(tmp_e.(fns{i}));
                        break;
                    end
                end
            end
        end
        
        % 2. 读取坐标轴 (取自第一阶Terms或Exact)
        tmp_coord = load(TermsFile, 'lon', 'lat');
        lon = wrapTo180(tmp_coord.lon);
        lat = tmp_coord.lat;
        [lon, sort_idx] = sort(lon); % 将 0~360 排序为 -180~180
        
        Exact_TSLA = Exact_TSLA(sort_idx, :, :); % 对齐经度
        [Nx, Ny, Nt] = size(Exact_TSLA);
        
        % 3. 读取并累加 Terms 分阶数据 (Lon x Lat x Time x Order)
        % 注意 Python 显示在后面，Matlab 读出来是 (Lon, Lat, Time, Order)
        tmp_t = load(TermsFile, 'TSLA_Terms_AllOrders');
        Terms_4D = single(tmp_t.TSLA_Terms_AllOrders(sort_idx, :, :, :));
        
        % 构建有效值掩码：Exact不为NaN，且该阶Term不为NaN
        mask_valid = ~isnan(Exact_TSLA) & ~all(isnan(Terms_4D), 4) & (Exact_TSLA ~= 0);
        
        % 开始逐阶计算与绘图
        Cumulative_Term = zeros(Nx, Ny, Nt, 'single');
        
        for p_idx = 1:8
            term_n = single(Terms_4D(:,:,:,p_idx));
            term_n(isnan(term_n)) = 0; % 防止 NaN 传播
            Cumulative_Term = Cumulative_Term + term_n;
            
            % 如果是我们感兴趣的阶数 (2~8)
            plot_row = find(Orders2Plot == p_idx);
            if ~isempty(plot_row)
                % 计算残差 = Exact - Taylor_Sum_N
                Residual = Exact_TSLA - Cumulative_Term;
                Residual(~mask_valid) = NaN;
                
                % 计算时间 RMSE
                RMSE_Map = sqrt(squeeze(mean(Residual .^ 2, 3, 'omitnan')));
                
                % 定位 subplot: 7x3 网格 (行=阶数, 列=模型)
                sub_id = (plot_row - 1) * 3 + m_idx;
                subplot(N_Orders, 3, sub_id);
                
                % 画图
                m_proj('robinson', 'long', [-180 180], 'lat', [-85 85]);
                m_pcolor(lon, lat, double(RMSE_Map)');
                shading flat;
                m_coast('color', [0.3 0.3 0.3], 'linewidth', 1);
                
                caxis(caxis_range);
                colormap(cmap);
                
                % 标题和标签
                if plot_row == 1 % 第一行带上模型名字
                    title(sprintf('%s (Order 1~%d)', model, p_idx), 'FontSize', 12, 'FontWeight', 'bold');
                else
                    title(sprintf('1~%d RMS Error', p_idx), 'FontSize', 10);
                end
                
                if m_idx == 1 % 第一列带上Y轴标签
                    ylabel(sprintf('Trunc: %d', p_idx), 'FontWeight', 'bold', 'Visible', 'on');
                end
                
                % 取消边框显示
                set(gca, 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
            end
        end
    end
    
    % 添加总体 Colorbar
    hp4 = get(subplot(N_Orders, 3, N_Orders * 3 - 1), 'Position');
    colorbar('Position', [hp4(1)+hp4(3)*0.1  0.03  hp4(3)*0.8  0.015], ...
             'Orientation', 'horizontal', 'FontSize', 12);
    annotation('textbox', [0 0.03 1 0.05], 'String', 'RMSE of Residual (mm)', ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % 总体标题
    sgtitle(sprintf('Truncation Error (RMSE) against Exact TSLA - [%s State]', state), 'FontSize', 18, 'FontWeight', 'bold');
    
    % 保存图片
    SaveFigurePath = fullfile(OutputFigureDir, sprintf('Residual_Comparison_%s.png', state));
    exportgraphics(fig, SaveFigurePath, 'Resolution', 300);
    fprintf('===图片 %s 已保存===\n', SaveFigurePath);
    
    close(fig);
end

fprintf('==== 全部绘制任务完成 ====\n');
warning(ws);
