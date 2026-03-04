clc
clear
close all

%========================================05-15
% load('Argo_all_05_24.mat')
%2000m以上===================================
load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\Argo_2000a_05_24.mat')
EN4_c13 = EN4_c13_2000a;
IAP = IAP_2000a;
Ishii = Ishii_2000a;

% load('Argo_2000a_05_15.mat')
% IAP =IAP_2000a;
% EN4_c13 =EN4_c13_2000a;
% Ishii = Ishii_2000a; 

% load('Argo_2000b_05_15.mat')
% IAP =IAP_2000b;
% EN4_c13 = EN4_c13_2000b; 
% Ishii = Ishii_2000b; 

start_year = 2005;
end_year = 2024;
index_05_24 = find(time_argo>start_year & time_argo<end_year+1);
%重力近似================================================================================


% %========================================05-24
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\IAP\IAP_05_24.mat')
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\EN4_c13\EN4_c13_05_24.mat')
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\Ishii\Ishii_05_24.mat')
% load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\IAP\IAP_all_new.mat', 'time_argo')
% 
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\IAP\IAP_2000a_05_24.mat')
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\EN4_c13\EN4_c13_2000a_05_24.mat')
% % load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\Ishii\Ishii_2000a_05_24.mat')
% % IAP =IAP_2000a;
% % EN4_c13 =EN4_c13_2000a;
% % Ishii = Ishii_2000a; 
% 
% load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\IAP\IAP_2000b_05_24.mat')
% load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\EN4_c13\EN4_c13_2000b_05_24.mat')
% load('F:\PRO_matlab\PRO_ARGO\PRO_TEST_SSLA\Ishii\Ishii_2000b_05_24.mat')
% IAP =IAP_2000b;
% EN4_c13 = EN4_c13_2000b; 
% Ishii = Ishii_2000b; index_05_24
% 
% start_year = 2005;
% end_year = 2024;
% index_05_24 = find(time_argo>start_year & time_argo<end_year+1);


%重力近似================================================================================
clearvars -except start_year end_year index_05_24 time_argo IAP EN4_c13 Ishii


for k= 1:25
    grid_IAP = IAP(k).grids;
    grid2_IAP = permute(grid_IAP(:,:,index_05_24),[2 1 3])*1000;
    [ ~, ~, ~,~, ~, ~, ~, ~, Trend, ~,~, ~, ~,~, ~, ~] = gmt_harmonic_new(time_argo(index_05_24),[],grid2_IAP);
    grids_IAP(:,k) = Trend(:);

    grid_EN4_c13 = EN4_c13(k).grids;
    grid2_EN4_c13 = permute(grid_EN4_c13(:,:,index_05_24),[2 1 3])*1000;
    [ ~, ~, ~,~, ~, ~, ~, ~, Trend, ~,~, ~, ~,~, ~, ~] = gmt_harmonic_new(time_argo(index_05_24),[],grid2_EN4_c13);
    grids_EN4_c13(:,k) = Trend(:);

    grid_Ishii = Ishii(k).grids;
    grid2_Ishii = permute(grid_Ishii(:,:,index_05_24),[2 1 3])*1000;
    [ ~, ~, ~,~, ~, ~, ~, ~, Trend, ~,~, ~, ~,~, ~, ~] = gmt_harmonic_new(time_argo(index_05_24),[],grid2_Ishii);
    grids_Ishii(:,k) = Trend(:);
end

index_nan = isnan(grids_IAP(:,2));
grids_IAP(index_nan,1) = NaN;
clean_grids_IAP = rmmissing(grids_IAP);

index_nan = isnan(grids_EN4_c13(:,2));
grids_EN4_c13(index_nan,1) = NaN;
clean_grids_EN4_c13 = rmmissing(grids_EN4_c13);

index_nan = isnan(grids_Ishii(:,2));
grids_Ishii(index_nan,1) = NaN;
clean_grids_Ishii = rmmissing(grids_Ishii);

clearvars -except start_year end_year index_05_24 time_argo clean_grids_IAP clean_grids_EN4_c13 clean_grids_Ishii




% 数据准备（请根据您的实际数据修改这部分）
% clean_grids_EN4_c13 = your_EN4_data;  % N×25 矩阵
% clean_grids_IAP = your_IAP_data;      % N×25 矩阵  
% clean_grids_Ishii = your_Ishii_data;  % N×25 矩阵

% ========== 可自定义参数 ==========
% 离群值抽稀比例（可自由调节）
outlier_sampling_ratio = 0.2; % 30%的离群值显示，可根据需要调整

% 宽度设置（可自由调节）
outlier_width = 0.2;          % 离群值的水平分布宽度
box_width = 0.4;              % 箱体宽度
whisker_line_width = 1;       % 须线竖线宽度
percentile_line_width = 3;    % 分位线横线宽度
median_line_width = 2;        % 中位线宽度
outlier_marker_size = 15;     % 离群值圆圈大小

% 箱线图位置设置
box_spacing = 0.8;              % 箱线图之间的间距
box_centers_3 = [1, 1+box_spacing, 1+2*box_spacing];    % 三个箱线图的中心位置
box_centers_4 = [1, 1+box_spacing, 1+2*box_spacing, 1+3*box_spacing]; % 四个箱线图的中心位置

% 字体设置
xlabel_fontsize = 16;         % 横轴字体大小
ylabel_fontsize = 16;         % 纵轴字体大小
title_fontsize = 16;          % 标题字体大小
tick_fontsize = 14;           % 刻度标签字体大小
legend_fontsize = 14;         % 图例字体大小

% 颜色设置 - 顶刊质量配色
colors = [0.6 0 0; 0 0 0.6; 0 0.6 0]; % 深色调：EN4-深红，IAP-深蓝，Ishii-深绿 - 用于离群值
box_color = [0.8 0.9 1]; % 箱体颜色 - 浅蓝色背景
border_color = [0.3 0.3 0.3]; % 边框颜色 - 深灰色
median_color = [0.8 0 0]; % 中位线颜色 - 红色
mean_color = [0 0 0]; % 均值点颜色 - 黑色
whisker_color = [0.3 0.3 0.3]; % 须线颜色 - 深灰色
percentile_color = [0.3 0.3 0.3]; % 分位线颜色 - 深灰色

categories_3 = {'SSLA', 'TSLA', 'HSLA'};
categories_4 = {'TSLA(1)', 'TSLA(2)', 'HSLA(1)', 'HSLA(2)'};

% 提取前三列并剔除0.5%极端值（两端各0.25%）
remove_extremes = @(data) data(data >= prctile(data, 0.25) & data <= prctile(data, 99.75));


% 创建图形，包含4个子图
figure('Position', [100, 100, 1600, 900], 'Color', 'white');

% ==================== 子图1: 第1-3列 ====================
subplot(1, 4, 1);
hold on;

% 处理EN4数据 - 第1-3列
steric_EN4_1 = remove_extremes(clean_grids_EN4_c13(:, 1));
thermo_EN4_1 = remove_extremes(clean_grids_EN4_c13(:, 2));
halo_EN4_1 = remove_extremes(clean_grids_EN4_c13(:, 3));

% 处理IAP数据 - 第1-3列
steric_IAP_1 = remove_extremes(clean_grids_IAP(:, 1));
thermo_IAP_1 = remove_extremes(clean_grids_IAP(:, 2));
halo_IAP_1 = remove_extremes(clean_grids_IAP(:, 3));

% 处理Ishii数据 - 第1-3列
steric_Ishii_1 = remove_extremes(clean_grids_Ishii(:, 1));
thermo_Ishii_1 = remove_extremes(clean_grids_Ishii(:, 2));
halo_Ishii_1 = remove_extremes(clean_grids_Ishii(:, 3));

% 将三套数据的比容、热容、盐容分别合并
steric_combined_1 = [steric_EN4_1; steric_IAP_1; steric_Ishii_1];
thermo_combined_1 = [thermo_EN4_1; thermo_IAP_1; thermo_Ishii_1];
halo_combined_1 = [halo_EN4_1; halo_IAP_1; halo_Ishii_1];

% 基于融合数据计算2.5%-97.5%分位线
steric_lower_95_1 = prctile(steric_combined_1, 2.5);
steric_upper_95_1 = prctile(steric_combined_1, 97.5);
thermo_lower_95_1 = prctile(thermo_combined_1, 2.5);
thermo_upper_95_1 = prctile(thermo_combined_1, 97.5);
halo_lower_95_1 = prctile(halo_combined_1, 2.5);
halo_upper_95_1 = prctile(halo_combined_1, 97.5);

% 准备数据
all_data_1 = {steric_combined_1, thermo_combined_1, halo_combined_1};
all_lower_95_1 = [steric_lower_95_1, thermo_lower_95_1, halo_lower_95_1];
all_upper_95_1 = [steric_upper_95_1, thermo_upper_95_1, halo_upper_95_1];

% 第一步：先绘制离群值（在最底层）
steric_data_sets_1 = {steric_EN4_1, steric_IAP_1, steric_Ishii_1};
for dataset_idx = 1:3
    data = steric_data_sets_1{dataset_idx};
    outliers = data(data < all_lower_95_1(1) | data > all_upper_95_1(1));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(1) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

thermo_data_sets_1 = {thermo_EN4_1, thermo_IAP_1, thermo_Ishii_1};
for dataset_idx = 1:3
    data = thermo_data_sets_1{dataset_idx};
    outliers = data(data < all_lower_95_1(2) | data > all_upper_95_1(2));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(2) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

halo_data_sets_1 = {halo_EN4_1, halo_IAP_1, halo_Ishii_1};
for dataset_idx = 1:3
    data = halo_data_sets_1{dataset_idx};
    outliers = data(data < all_lower_95_1(3) | data > all_upper_95_1(3));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(3) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 第二步：绘制箱线图主体（在离群值之上）
for i = 1:3
    data = all_data_1{i};
    center = box_centers_3(i);
    
    q1 = prctile(data, 25);
    q3 = prctile(data, 75);
    median_val = median(data);
    mean_val = mean(data);
    
    x_left = center - box_width/2;
    x_right = center + box_width/2;
    
    patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
    plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);
    
    lower_whisker = all_lower_95_1(i);
    upper_whisker = all_upper_95_1(i);
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
    plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
    scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% 第三步：最后绘制分位线横线（在最上层，盖住所有元素）
for i = 1:3
    center = box_centers_3(i);
    lower_whisker = all_lower_95_1(i);
    upper_whisker = all_upper_95_1(i);
    data = all_data_1{i};
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
    plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);
end

ylabel('Trend bais (mm/yr)', 'FontSize', ylabel_fontsize, 'FontWeight', 'bold');
set(gca, 'XTick', box_centers_3, 'XTickLabel', categories_3, 'FontSize', tick_fontsize, 'FontWeight', 'bold');
grid on;
grid minor;

% 添加子图标签 a)
text(0.02, 0.985, '(13) Gravitational approximation', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold',  ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% ==================== 子图2: 第4-12列 ====================
subplot(1, 4, 2);
hold on;

% 处理EN4数据 - 第4-12列（4-6比容，7-9热容，10-12盐容）
steric_EN4_2 = remove_extremes(clean_grids_EN4_c13(:, 5));
thermo_EN4_2 = remove_extremes(clean_grids_EN4_c13(:, 8));
halo_EN4_2 = remove_extremes(clean_grids_EN4_c13(:, 11));

% 处理IAP数据 - 第4-12列
steric_IAP_2 = remove_extremes(clean_grids_IAP(:, 5));
thermo_IAP_2 = remove_extremes(clean_grids_IAP(:, 8));
halo_IAP_2 = remove_extremes(clean_grids_IAP(:, 11));

% 处理Ishii数据 - 第4-12列
steric_Ishii_2 = remove_extremes(clean_grids_Ishii(:, 5));
thermo_Ishii_2 = remove_extremes(clean_grids_Ishii(:, 8));
halo_Ishii_2 = remove_extremes(clean_grids_Ishii(:, 11));

% 将三套数据的比容、热容、盐容分别合并
steric_combined_2 = [steric_EN4_2; steric_IAP_2; steric_Ishii_2];
thermo_combined_2 = [thermo_EN4_2; thermo_IAP_2; thermo_Ishii_2];
halo_combined_2 = [halo_EN4_2; halo_IAP_2; halo_Ishii_2];

% 基于融合数据计算2.5%-97.5%分位线
steric_lower_95_2 = prctile(steric_combined_2, 2.5);
steric_upper_95_2 = prctile(steric_combined_2, 97.5);
thermo_lower_95_2 = prctile(thermo_combined_2, 2.5);
thermo_upper_95_2 = prctile(thermo_combined_2, 97.5);
halo_lower_95_2 = prctile(halo_combined_2, 2.5);
halo_upper_95_2 = prctile(halo_combined_2, 97.5);

% 准备数据
all_data_2 = {steric_combined_2, thermo_combined_2, halo_combined_2};
all_lower_95_2 = [steric_lower_95_2, thermo_lower_95_2, halo_lower_95_2];
all_upper_95_2 = [steric_upper_95_2, thermo_upper_95_2, halo_upper_95_2];

% 第一步：先绘制离群值（在最底层）
steric_data_sets_2 = {steric_EN4_2, steric_IAP_2, steric_Ishii_2};
for dataset_idx = 1:3
    data = steric_data_sets_2{dataset_idx};
    outliers = data(data < all_lower_95_2(1) | data > all_upper_95_2(1));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(1) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

thermo_data_sets_2 = {thermo_EN4_2, thermo_IAP_2, thermo_Ishii_2};
for dataset_idx = 1:3
    data = thermo_data_sets_2{dataset_idx};
    outliers = data(data < all_lower_95_2(2) | data > all_upper_95_2(2));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(2) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

halo_data_sets_2 = {halo_EN4_2, halo_IAP_2, halo_Ishii_2};
for dataset_idx = 1:3
    data = halo_data_sets_2{dataset_idx};
    outliers = data(data < all_lower_95_2(3) | data > all_upper_95_2(3));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(3) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 第二步：绘制箱线图主体（在离群值之上）
for i = 1:3
    data = all_data_2{i};
    center = box_centers_3(i);
    
    q1 = prctile(data, 25);
    q3 = prctile(data, 75);
    median_val = median(data);
    mean_val = mean(data);
    
    x_left = center - box_width/2;
    x_right = center + box_width/2;
    
    patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
    plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);
    
    lower_whisker = all_lower_95_2(i);
    upper_whisker = all_upper_95_2(i);
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
    plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
    scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% 第三步：最后绘制分位线横线（在最上层，盖住所有元素）
for i = 1:3
    center = box_centers_3(i);
    lower_whisker = all_lower_95_2(i);
    upper_whisker = all_upper_95_2(i);
    data = all_data_2{i};
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
    plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);
end

set(gca, 'XTick', box_centers_3, 'XTickLabel', categories_3, 'FontSize', tick_fontsize, 'FontWeight', 'bold');
grid on;
grid minor;

% 添加子图标签 b)
text(0.02, 0.985, '(14) Density approximation', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% ==================== 子图3: 第13-15列（双y轴） ====================
subplot(1, 4, 3);
hold on;

% 处理EN4数据 - 第13-15列
steric_EN4_3 = remove_extremes(clean_grids_EN4_c13(:, 13));
thermo_EN4_3 = remove_extremes(clean_grids_EN4_c13(:, 14));
halo_EN4_3 = remove_extremes(clean_grids_EN4_c13(:, 15));

% 处理IAP数据 - 第13-15列
steric_IAP_3 = remove_extremes(clean_grids_IAP(:, 13));
thermo_IAP_3 = remove_extremes(clean_grids_IAP(:, 14));
halo_IAP_3 = remove_extremes(clean_grids_IAP(:, 15));

% 处理Ishii数据 - 第13-15列
steric_Ishii_3 = remove_extremes(clean_grids_Ishii(:, 13));
thermo_Ishii_3 = remove_extremes(clean_grids_Ishii(:, 14));
halo_Ishii_3 = remove_extremes(clean_grids_Ishii(:, 15));

% 将三套数据的比容、热容、盐容分别合并
steric_combined_3 = [steric_EN4_3; steric_IAP_3; steric_Ishii_3];
thermo_combined_3 = [thermo_EN4_3; thermo_IAP_3; thermo_Ishii_3];
halo_combined_3 = [halo_EN4_3; halo_IAP_3; halo_Ishii_3];

% 基于融合数据计算2.5%-97.5%分位线
steric_lower_95_3 = prctile(steric_combined_3, 2.5);
steric_upper_95_3 = prctile(steric_combined_3, 97.5);
thermo_lower_95_3 = prctile(thermo_combined_3, 2.5);
thermo_upper_95_3 = prctile(thermo_combined_3, 97.5);
halo_lower_95_3 = prctile(halo_combined_3, 2.5);
halo_upper_95_3 = prctile(halo_combined_3, 97.5);

% 准备数据
all_data_3 = {steric_combined_3, thermo_combined_3, halo_combined_3};
all_lower_95_3 = [steric_lower_95_3, thermo_lower_95_3, halo_lower_95_3];
all_upper_95_3 = [steric_upper_95_3, thermo_upper_95_3, halo_upper_95_3];

% 颜色设置 - 顶刊质量双轴配色
left_color = [0 0.4 0.8];   % 左侧y轴颜色 - 蓝色
right_color = [0.8 0.2 0];  % 右侧y轴颜色 - 红色

% 创建左侧y轴（比容和盐容）
yyaxis left;
hold on;

% 绘制比容（左侧y轴）
data = all_data_3{1};
center = box_centers_3(1);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_3(1);
upper_whisker = all_upper_95_3(1);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制比容的离群值（左侧y轴）
steric_data_sets_3 = {steric_EN4_3, steric_IAP_3, steric_Ishii_3};
for dataset_idx = 1:3
    data = steric_data_sets_3{dataset_idx};
    outliers = data(data < all_lower_95_3(1) | data > all_upper_95_3(1));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(1) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 绘制盐容（左侧y轴）
data = all_data_3{3};
center = box_centers_3(3);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_3(3);
upper_whisker = all_upper_95_3(3);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制盐容的离群值（左侧y轴）
halo_data_sets_3 = {halo_EN4_3, halo_IAP_3, halo_Ishii_3};
for dataset_idx = 1:3
    data = halo_data_sets_3{dataset_idx};
    outliers = data(data < all_lower_95_3(3) | data > all_upper_95_3(3));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(3) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 设置左侧y轴
set(gca, 'YColor', left_color, 'FontWeight', 'bold');

% 创建右侧y轴（热容）
yyaxis right;
hold on;

% 绘制热容（右侧y轴）
data = all_data_3{2};
center = box_centers_3(2);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_3(2);
upper_whisker = all_upper_95_3(2);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制热容的离群值（右侧y轴）
thermo_data_sets_3 = {thermo_EN4_3, thermo_IAP_3, thermo_Ishii_3};
for dataset_idx = 1:3
    data = thermo_data_sets_3{dataset_idx};
    outliers = data(data < all_lower_95_3(2) | data > all_upper_95_3(2));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_3(2) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 设置右侧y轴
set(gca, 'YColor', right_color, 'FontWeight', 'bold');

% 最后绘制分位线横线（在最上层，盖住所有元素）
yyaxis left;
for i = [1, 3]  % 比容和盐容
    center = box_centers_3(i);
    lower_whisker = all_lower_95_3(i);
    upper_whisker = all_upper_95_3(i);
    data = all_data_3{i};
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
    plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);
end

yyaxis right;
center = box_centers_3(2);  % 热容
lower_whisker = all_lower_95_3(2);
upper_whisker = all_upper_95_3(2);
data = all_data_3{2};
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);

% 设置横轴标签颜色
ax = gca;
ax.XTick = box_centers_3;
ax.XTickLabel = categories_3;
ax.XTickLabel{1} = ['\color[rgb]{' num2str(left_color) '}' categories_3{1}];
ax.XTickLabel{2} = ['\color[rgb]{' num2str(right_color) '}' categories_3{2}];
ax.XTickLabel{3} = ['\color[rgb]{' num2str(left_color) '}' categories_3{3}];
set(gca, 'FontSize', tick_fontsize, 'FontWeight', 'bold');
grid on;
grid minor;

% 添加子图标签 c)
text(0.02, 0.985, '(15) Reference selection', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% ==================== 子图4: 第16-19列（双y轴） ====================

subplot(1, 4, 4);
hold on;

% 处理EN4数据 - 第16-19列
thermo1_EN4_4 = remove_extremes(clean_grids_EN4_c13(:, 16));
thermo2_EN4_4 = remove_extremes(clean_grids_EN4_c13(:, 17));
halo1_EN4_4 = remove_extremes(clean_grids_EN4_c13(:, 18));
halo2_EN4_4 = remove_extremes(clean_grids_EN4_c13(:, 19));

% 处理IAP数据 - 第16-19列
thermo1_IAP_4 = remove_extremes(clean_grids_IAP(:, 16));
thermo2_IAP_4 = remove_extremes(clean_grids_IAP(:, 17));
halo1_IAP_4 = remove_extremes(clean_grids_IAP(:, 18));
halo2_IAP_4 = remove_extremes(clean_grids_IAP(:, 19));

% 处理Ishii数据 - 第16-19列
thermo1_Ishii_4 = remove_extremes(clean_grids_Ishii(:, 16));
thermo2_Ishii_4 = remove_extremes(clean_grids_Ishii(:, 17));
halo1_Ishii_4 = remove_extremes(clean_grids_Ishii(:, 18));
halo2_Ishii_4 = remove_extremes(clean_grids_Ishii(:, 19));

% 将三套数据分别合并
thermo1_combined_4 = [thermo1_EN4_4; thermo1_IAP_4; thermo1_Ishii_4];
thermo2_combined_4 = [thermo2_EN4_4; thermo2_IAP_4; thermo2_Ishii_4];
halo1_combined_4 = [halo1_EN4_4; halo1_IAP_4; halo1_Ishii_4];
halo2_combined_4 = [halo2_EN4_4; halo2_IAP_4; halo2_Ishii_4];

% 基于融合数据计算2.5%-97.5%分位线
thermo1_lower_95_4 = prctile(thermo1_combined_4, 2.5);
thermo1_upper_95_4 = prctile(thermo1_combined_4, 97.5);
thermo2_lower_95_4 = prctile(thermo2_combined_4, 2.5);
thermo2_upper_95_4 = prctile(thermo2_combined_4, 97.5);
halo1_lower_95_4 = prctile(halo1_combined_4, 2.5);
halo1_upper_95_4 = prctile(halo1_combined_4, 97.5);
halo2_lower_95_4 = prctile(halo2_combined_4, 2.5);
halo2_upper_95_4 = prctile(halo2_combined_4, 97.5);

% 准备数据
all_data_4 = {thermo1_combined_4, thermo2_combined_4, halo1_combined_4, halo2_combined_4};
all_lower_95_4 = [thermo1_lower_95_4, thermo2_lower_95_4, halo1_lower_95_4, halo2_lower_95_4];
all_upper_95_4 = [thermo1_upper_95_4, thermo2_upper_95_4, halo1_upper_95_4, halo2_upper_95_4];

% 颜色设置 - 顶刊质量双轴配色
thermo_color = [0.8 0.2 0];  % 热容颜色 - 红色
halo_color = [0 0.4 0.8];    % 盐容颜色 - 蓝色

% 创建左侧y轴（热容一阶）
yyaxis left;
hold on;

% 绘制热容一阶（左侧y轴）
data = all_data_4{1};
center = box_centers_4(1);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_4(1);
upper_whisker = all_upper_95_4(1);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制热容一阶的离群值
thermo1_data_sets_4 = {thermo1_EN4_4, thermo1_IAP_4, thermo1_Ishii_4};
for dataset_idx = 1:3
    data = thermo1_data_sets_4{dataset_idx};
    outliers = data(data < all_lower_95_4(1) | data > all_upper_95_4(1));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_4(1) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 设置左侧y轴范围 - 基于热容一阶数据
data_range_left = [min(thermo1_combined_4), max(thermo1_combined_4)];
padding_left = (data_range_left(2) - data_range_left(1)) * 0.1; % 10%的边距
ylim_left = [data_range_left(1)-padding_left, data_range_left(2)+padding_left];
ylim(ylim_left);
set(gca, 'YColor', thermo_color, 'FontWeight', 'bold');

% 创建右侧y轴（热容二阶、盐容一阶、盐容二阶）
yyaxis right;
hold on;

% 绘制热容二阶（右侧y轴）
data = all_data_4{2};
center = box_centers_4(2);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_4(2);
upper_whisker = all_upper_95_4(2);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制热容二阶的离群值
thermo2_data_sets_4 = {thermo2_EN4_4, thermo2_IAP_4, thermo2_Ishii_4};
for dataset_idx = 1:3
    data = thermo2_data_sets_4{dataset_idx};
    outliers = data(data < all_lower_95_4(2) | data > all_upper_95_4(2));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_4(2) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 绘制盐容一阶（右侧y轴）
data = all_data_4{3};
center = box_centers_4(3);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_4(3);
upper_whisker = all_upper_95_4(3);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制盐容一阶的离群值
halo1_data_sets_4 = {halo1_EN4_4, halo1_IAP_4, halo1_Ishii_4};
for dataset_idx = 1:3
    data = halo1_data_sets_4{dataset_idx};
    outliers = data(data < all_lower_95_4(3) | data > all_upper_95_4(3));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_4(3) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 绘制盐容二阶（右侧y轴）
data = all_data_4{4};
center = box_centers_4(4);

q1 = prctile(data, 25);
q3 = prctile(data, 75);
median_val = median(data);
mean_val = mean(data);

x_left = center - box_width/2;
x_right = center + box_width/2;

patch([x_left, x_right, x_right, x_left], [q1, q1, q3, q3], box_color, 'FaceAlpha', 0.7, 'EdgeColor', border_color, 'LineWidth', 1.5);
plot([x_left, x_right], [median_val, median_val], 'Color', median_color, 'LineWidth', median_line_width);

lower_whisker = all_lower_95_4(4);
upper_whisker = all_upper_95_4(4);
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center, center], [q1, whisker_min], 'Color', whisker_color, 'LineWidth', whisker_line_width);
plot([center, center], [q3, whisker_max], 'Color', whisker_color, 'LineWidth', whisker_line_width);
scatter(center, mean_val, 100, mean_color, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 绘制盐容二阶的离群值
halo2_data_sets_4 = {halo2_EN4_4, halo2_IAP_4, halo2_Ishii_4};
for dataset_idx = 1:3
    data = halo2_data_sets_4{dataset_idx};
    outliers = data(data < all_lower_95_4(4) | data > all_upper_95_4(4));
    if ~isempty(outliers)
        n_outliers = length(outliers);
        n_display = max(1, round(n_outliers * outlier_sampling_ratio));
        if n_display < n_outliers
            idx = randperm(n_outliers, n_display);
            outliers = outliers(idx);
        end
        x_pos = box_centers_4(4) + (rand(size(outliers)) - 0.5) * outlier_width;
        scatter(x_pos, outliers, outlier_marker_size, colors(dataset_idx, :), 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

% 设置右侧y轴范围 - 基于热容二阶、盐容一阶、盐容二阶数据的整体范围
all_right_data = [thermo2_combined_4; halo1_combined_4; halo2_combined_4];
data_range_right = [min(all_right_data), max(all_right_data)];
padding_right = (data_range_right(2) - data_range_right(1)) * 0.1; % 10%的边距
ylim_right = [data_range_right(1)-padding_right, data_range_right(2)+padding_right];
ylim(ylim_right);
set(gca, 'YColor', halo_color, 'FontWeight', 'bold');

% 最后绘制分位线横线（在最上层，盖住所有元素）
yyaxis left;
center = box_centers_4(1);  % 热容一阶
lower_whisker = all_lower_95_4(1);
upper_whisker = all_upper_95_4(1);
data = all_data_4{1};
data_in_range = data(data >= lower_whisker & data <= upper_whisker);
whisker_min = min(data_in_range);
whisker_max = max(data_in_range);

plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);

yyaxis right;
for i = [2, 3, 4]  % 热容二阶、盐容一阶、盐容二阶
    center = box_centers_4(i);
    lower_whisker = all_lower_95_4(i);
    upper_whisker = all_upper_95_4(i);
    data = all_data_4{i};
    data_in_range = data(data >= lower_whisker & data <= upper_whisker);
    whisker_min = min(data_in_range);
    whisker_max = max(data_in_range);
    
    plot([center-0.1, center+0.1], [whisker_min, whisker_min], 'Color', percentile_color, 'LineWidth', percentile_line_width);
    plot([center-0.1, center+0.1], [whisker_max, whisker_max], 'Color', percentile_color, 'LineWidth', percentile_line_width);
end

% 设置横轴标签颜色
ax = gca;
ax.XTick = box_centers_4;
ax.XTickLabel = categories_4;
ax.XTickLabel{1} = ['\color[rgb]{' num2str(thermo_color) '}' categories_4{1}];
ax.XTickLabel{2} = ['\color[rgb]{' num2str(halo_color) '}' categories_4{2}];
ax.XTickLabel{3} = ['\color[rgb]{' num2str(halo_color) '}' categories_4{3}];
ax.XTickLabel{4} = ['\color[rgb]{' num2str(halo_color) '}' categories_4{4}];
set(gca, 'FontSize', tick_fontsize, 'FontWeight', 'bold');
grid on;
grid minor;

% 添加子图标签 d)
text(0.02, 0.985, '(16) Low-order approximation', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% 调整子图间距
set(gcf, 'Position', [100, 100, 1600, 900]); % 调整图形大小

% 在主图上添加统一图例（放在顶部）
legend_elements = [
    plot(NaN, NaN, 'Color', median_color, 'LineWidth', 2, 'DisplayName', 'Median line'), ...
    plot(NaN, NaN, 'ko', 'MarkerFaceColor', mean_color, 'MarkerSize', 8, 'DisplayName', 'Mean value'), ...
    plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', percentile_color, 'DisplayName', '2.5% ~ 97.5%'), ...
    plot(NaN, NaN, 'ro', 'MarkerFaceColor', colors(1, :), 'MarkerSize', 8, 'DisplayName', 'Outliers-EN4'), ...
    plot(NaN, NaN, 'bo', 'MarkerFaceColor', colors(2, :), 'MarkerSize', 8, 'DisplayName', 'Outliers-IAP'), ...
    plot(NaN, NaN, 'go', 'MarkerFaceColor', colors(3, :), 'MarkerSize', 8, 'DisplayName', 'Outliers-Ishii')
];

% 在主图上添加图例
% 精细调整
hL = legend(legend_elements, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
       'NumColumns', 6, 'FontSize', legend_fontsize, 'Box', 'off', 'FontWeight', 'bold');

% 让图形和文字几乎挨着
hL.ItemTokenSize = [50, 18];   % 非常小的值，让图形和文字紧贴

% 保持图例项之间有足够间距
if isprop(hL, 'ColumnSpacing')
    hL.ColumnSpacing = 100;    % 较大的列间距
end

% 给图例足够的宽度来容纳间距
hL.Position = [0.20, 0.95, 0.94, 0.05];

% 调整子图位置使其更紧凑
for i = 1:4
    h = subplot(1, 4, i);
    box on;
    pos = h.Position;
    if i == 1
        pos(1) = 0.05; % 第一个子图左移
    elseif i == 2
        pos(1) = 0.28; % 第二个子图左移
    elseif i == 3
        pos(1) = 0.51; % 第三个子图左移
    else
        pos(1) = 0.75; % 第四个子图左移
    end
    pos(3) = 0.18; % 统一子图宽度
    h.Position = pos;
end


%% ========== 输出箱线图统计结果（简洁版） ==========

% 创建输出文件夹
output_dir = 'boxplot_statistics';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 获取当前日期时间作为文件名的一部分
current_time = datestr(now, 'yyyy-mm-dd_HH-MM');
output_filename = fullfile(output_dir, ['boxplot_statistics_', current_time, '.txt']);

% 打开文件准备写入
fid = fopen(output_filename, 'w');

if fid == -1
    error('无法创建输出文件');
end

%% 子图1统计结果
fprintf(fid, '【子图1：第1-3列 - SSLA, TSLA, HSLA】\n');
fprintf(fid, '%-10s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    '类别', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');
fprintf(fid, '%s\n', repmat('-', 145, 1));

for i = 1:3
    % 融合数据（剔除极端值后）
    data = all_data_1{i};
    category = categories_3{i};
    
    % 获取原始数据（剔除极端值前）
    if i == 1 % SSLA
        orig_data = [clean_grids_EN4_c13(:,1); clean_grids_IAP(:,1); clean_grids_Ishii(:,1)];
    elseif i == 2 % TSLA
        orig_data = [clean_grids_EN4_c13(:,2); clean_grids_IAP(:,2); clean_grids_Ishii(:,2)];
    else % HSLA
        orig_data = [clean_grids_EN4_c13(:,3); clean_grids_IAP(:,3); clean_grids_Ishii(:,3)];
    end
    
    orig_min = min(orig_data);
    orig_max = max(orig_data);
    filtered_min = min(data);
    filtered_max = max(data);
    p2_5 = prctile(data, 2.5);
    p25 = prctile(data, 25);
    median_val = median(data);
    p75 = prctile(data, 75);
    p97_5 = prctile(data, 97.5);
    mean_val = mean(data);
    
    fprintf(fid, '%-10s %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        category, orig_min, orig_max, filtered_min, filtered_max, p2_5, p25, median_val, p75, p97_5, mean_val);
end

fprintf(fid, '\n\n');

%% 子图2统计结果
fprintf(fid, '【子图2：第4-12列（第5,8,11列）】\n');
fprintf(fid, '%-10s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    '类别', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');
fprintf(fid, '%s\n', repmat('-', 145, 1));

col_indices = [5, 8, 11]; % 第5,8,11列
for i = 1:3
    % 融合数据（剔除极端值后）
    data = all_data_2{i};
    category = categories_3{i};
    
    % 获取原始数据（剔除极端值前）
    col_idx = col_indices(i);
    orig_data = [clean_grids_EN4_c13(:,col_idx); clean_grids_IAP(:,col_idx); clean_grids_Ishii(:,col_idx)];
    
    orig_min = min(orig_data);
    orig_max = max(orig_data);
    filtered_min = min(data);
    filtered_max = max(data);
    p2_5 = prctile(data, 2.5);
    p25 = prctile(data, 25);
    median_val = median(data);
    p75 = prctile(data, 75);
    p97_5 = prctile(data, 97.5);
    mean_val = mean(data);
    
    fprintf(fid, '%-10s %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        category, orig_min, orig_max, filtered_min, filtered_max, p2_5, p25, median_val, p75, p97_5, mean_val);
end

fprintf(fid, '\n\n');

%% 子图3统计结果（双Y轴）
fprintf(fid, '【子图3：第13-15列（双Y轴）】\n');
fprintf(fid, '%-15s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    '类别(Y轴)', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');
fprintf(fid, '%s\n', repmat('-', 150, 1));

col_indices_3 = [13, 14, 15]; % 第13,14,15列
categories_3_with_axis = {'SSLA(左轴)', 'TSLA(右轴)', 'HSLA(左轴)'};

for i = 1:3
    % 融合数据（剔除极端值后）
    data = all_data_3{i};
    category = categories_3_with_axis{i};
    
    % 获取原始数据（剔除极端值前）
    col_idx = col_indices_3(i);
    orig_data = [clean_grids_EN4_c13(:,col_idx); clean_grids_IAP(:,col_idx); clean_grids_Ishii(:,col_idx)];
    
    orig_min = min(orig_data);
    orig_max = max(orig_data);
    filtered_min = min(data);
    filtered_max = max(data);
    p2_5 = prctile(data, 2.5);
    p25 = prctile(data, 25);
    median_val = median(data);
    p75 = prctile(data, 75);
    p97_5 = prctile(data, 97.5);
    mean_val = mean(data);
    
    fprintf(fid, '%-15s %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        category, orig_min, orig_max, filtered_min, filtered_max, p2_5, p25, median_val, p75, p97_5, mean_val);
end

fprintf(fid, '\n\n');

%% 子图4统计结果（双Y轴）
fprintf(fid, '【子图4：第16-19列】\n');
fprintf(fid, '%-15s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    '类别(Y轴)', '原始最小', '原始最大', '剔除后最小', '剔除后最大', '2.5%', '25%', '中位数', '75%', '97.5%', '均值');
fprintf(fid, '%s\n', repmat('-', 150, 1));

col_indices_4 = [16, 17, 18, 19]; % 第16-19列
categories_4_with_axis = {'TSLA(M1)(左轴)', 'TSLA(M2)(右轴)', 'HSLA(M1)(右轴)', 'HSLA(M2)(右轴)'};

for i = 1:4
    % 融合数据（剔除极端值后）
    data = all_data_4{i};
    category = categories_4_with_axis{i};
    
    % 获取原始数据（剔除极端值前）
    col_idx = col_indices_4(i);
    orig_data = [clean_grids_EN4_c13(:,col_idx); clean_grids_IAP(:,col_idx); clean_grids_Ishii(:,col_idx)];
    
    orig_min = min(orig_data);
    orig_max = max(orig_data);
    filtered_min = min(data);
    filtered_max = max(data);
    p2_5 = prctile(data, 2.5);
    p25 = prctile(data, 25);
    median_val = median(data);
    p75 = prctile(data, 75);
    p97_5 = prctile(data, 97.5);
    mean_val = mean(data);
    
    fprintf(fid, '%-15s %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        category, orig_min, orig_max, filtered_min, filtered_max, p2_5, p25, median_val, p75, p97_5, mean_val);
end

%% 关闭文件
fclose(fid);

%% 在命令行显示完成信息
fprintf('统计结果已保存到: %s\n', output_filename);

%% 保持图形打开
uiwait;