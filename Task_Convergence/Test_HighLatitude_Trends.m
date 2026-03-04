function Test_HighLatitude_Trends()
% 功能: 分析高纬度地区的二阶及以上趋势
% 输入: 无
% 输出: 高纬度地区趋势统计信息和可视化结果

clc; close all;
fprintf('\n======================================================\n');
fprintf('🧪 启动高纬度趋势分析测试\n');
fprintf('======================================================\n');

%% 1. 数据加载
TrendDir = 'd:\\work\\Task_Convergence\\Trend_Results';
datasets = {'EN4', 'IAP', 'Ishii'};
states = {'Average', 'StdRef'};

%% 2. 高纬度定义 (纬度绝对值大于60度)
high_lat_threshold = 75;

%% 3. 分析每个数据集和状态
for d = 1:length(datasets)
    dataset_name = datasets{d};
    fprintf('\n>> 分析 %s 数据...\n', dataset_name);
    
    for s = 1:length(states)
        state_name = states{s};
        fprintf('   - 状态: %s\n', state_name);
        
        % 加载趋势文件
        file_name = fullfile(TrendDir, sprintf('%s_%s_Trends.mat', dataset_name, state_name));
        if ~exist(file_name, 'file')
            fprintf('   ⚠️ 文件不存在: %s\n', file_name);
            continue;
        end
        
        data = load(file_name);
        
        % 提取纬度和趋势数据
        lat = data.Lat;
        lon = data.Lon;
        
        % 计算高纬度掩码
        high_lat_mask = abs(lat) > high_lat_threshold;
        n_high_lat = sum(high_lat_mask);
        fprintf('   - 高纬度网格点数量: %d\n', n_high_lat);
        
        if n_high_lat == 0
            fprintf('   ⚠️ 无高纬度数据\n');
            continue;
        end
        
        % 提取趋势数据
        trend_TSLA = data.trend_TSLA;  % T项趋势
        trend_SSLA = data.trend_SSLA;  % S项趋势
        trend_HSLA = data.trend_HSLA;  % 混合项趋势
        
        % 计算各阶趋势的高纬度平均值
        fprintf('   - 高纬度趋势统计 (单位: mm/yr):\n');
        
        % T项趋势 (T1-T3)
        for i = 1:3
            t_trend = squeeze(trend_TSLA(:, :, i));
            t_high_lat = t_trend(high_lat_mask, :);
            t_mean = nanmean(t_high_lat(:));
            t_std = nanstd(t_high_lat(:));
            fprintf('     T%d: 均值 = %.4f, 标准差 = %.4f\n', i, t_mean, t_std);
        end
        
        % S项趋势 (S1-S3)
        for i = 1:3
            s_trend = squeeze(trend_SSLA(:, :, i));
            s_high_lat = s_trend(high_lat_mask, :);
            s_mean = nanmean(s_high_lat(:));
            s_std = nanstd(s_high_lat(:));
            fprintf('     S%d: 均值 = %.4f, 标准差 = %.4f\n', i, s_mean, s_std);
        end
        
        % 混合项趋势 (H1-H8)
        for i = 1:8
            h_trend = squeeze(trend_HSLA(:, :, i));
            h_high_lat = h_trend(high_lat_mask, :);
            h_mean = nanmean(h_high_lat(:));
            h_std = nanstd(h_high_lat(:));
            fprintf('     H%d: 均值 = %.4f, 标准差 = %.4f\n', i, h_mean, h_std);
        end
        
        % 计算二阶及以上项的总趋势
        t_second_order = sum(trend_TSLA(:, :, 2:3), 3);
        s_second_order = sum(trend_SSLA(:, :, 2:3), 3);
        h_second_order = sum(trend_HSLA(:, :, 2:8), 3);  % 混合项从H2开始算二阶及以上
        
        t_second_high = t_second_order(high_lat_mask, :);
        s_second_high = s_second_order(high_lat_mask, :);
        h_second_high = h_second_order(high_lat_mask, :);
        
        total_second_order = t_second_high + s_second_high + h_second_high;
        
        fprintf('   - 二阶及以上总趋势: 均值 = %.4f, 标准差 = %.4f\n', ...
            nanmean(total_second_order(:)), nanstd(total_second_order(:)));
        
        % 可视化高纬度趋势
        figure('Position', [50, 50, 1200, 800], 'Color', 'w', 'Name', sprintf('%s %s 高纬度趋势', dataset_name, state_name));
        
        % 绘制T项趋势
        subplot(3, 1, 1);
        t_total = sum(trend_TSLA, 3);
        t_total_high = t_total(high_lat_mask, :);
        hist(t_total_high(:), 50);
        title('T项总趋势 (高纬度地区)', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('趋势 (mm/yr)', 'FontSize', 12);
        ylabel('频率', 'FontSize', 12);
        grid on;
        
        % 绘制S项趋势
        subplot(3, 1, 2);
        s_total = sum(trend_SSLA, 3);
        s_total_high = s_total(high_lat_mask, :);
        hist(s_total_high(:), 50);
        title('S项总趋势 (高纬度地区)', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('趋势 (mm/yr)', 'FontSize', 12);
        ylabel('频率', 'FontSize', 12);
        grid on;
        
        % 绘制混合项趋势
        subplot(3, 1, 3);
        h_total = sum(trend_HSLA, 3);
        h_total_high = h_total(high_lat_mask, :);
        hist(h_total_high(:), 50);
        title('混合项总趋势 (高纬度地区)', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('趋势 (mm/yr)', 'FontSize', 12);
        ylabel('频率', 'FontSize', 12);
        grid on;
        
        tightfig;
        
        % 保存图表
        saveas(gcf, sprintf('d:\\work\\Task_Convergence\\HighLatitude_Trends_%s_%s.png', dataset_name, state_name));
        fprintf('   - 图表已保存\n');
    end
end

fprintf('\n======================================================\n');
fprintf('✅ 高纬度趋势分析测试完成！\n');
fprintf('======================================================\n');
end

% 辅助函数: 调整图表布局
function tightfig()
    set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
    drawnow;
end