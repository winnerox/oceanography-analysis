function Calculate_Trends()
% 用最小二乘法一元线性回归 (OLS Linear Regression) 方法
% 计算三套数据的平均态和标准态的trend文件

fprintf('========== 计算泰勒展开项趋势 ==========\n\n');

% 基础目录
base_dir = 'D:\work';
output_dir = fullfile(base_dir, 'Task_Convergence', 'Trend_Results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 三套数据目录
data_dirs = {
    struct('name', 'EN4', 'path', fullfile(base_dir, 'EN4_TSLA_Terms'));
    struct('name', 'IAP', 'path', fullfile(base_dir, 'IAP_TSLA_Terms'));
    struct('name', 'Ishii', 'path', fullfile(base_dir, 'Ishii_mat_data'));
};

% 两种状态
states = {'Average', 'StdRef'};
state_cn = {'平均态', '标准态'};

% 时间向量（假设2005-2024年，每月一个数据点，共240个月）
years = 2005:2024;
months = 1:12;
time_years = zeros(length(years) * length(months), 1);
for i = 1:length(years)
    for j = 1:length(months)
        idx = (i-1)*12 + j;
        time_years(idx) = years(i) + (months(j)-1)/12;
    end
end

% 主循环：计算每套数据、每种状态、每阶的趋势
for d_idx = 1:length(data_dirs)
    dataset = data_dirs{d_idx};
    fprintf('\n处理数据集: %s\n', dataset.name);
    
    for s_idx = 1:length(states)
        state = states{s_idx};
        fprintf('  状态: %s (%s)\n', state_cn{s_idx}, state);
        
        % 构建文件名
        tsla_file = fullfile(dataset.path, sprintf('%s_TSLA_Terms_1to8_%s.mat', dataset.name, state));
        s_file = fullfile(dataset.path, sprintf('%s_S_Terms_1to8_%s.mat', dataset.name, state));
        cross_file = fullfile(dataset.path, sprintf('%s_Cross_Terms_1to8_%s.mat', dataset.name, state));
        
        % 检查文件是否存在
        if ~exist(tsla_file, 'file')
            fprintf('    警告: 文件不存在 %s\n', tsla_file);
            continue;
        end
        
        % 加载数据
        fprintf('    加载数据...\n');
        try
            % 使用load函数加载MAT文件
            tsla_data = load(tsla_file);
            s_data = load(s_file);
            cross_data = load(cross_file);
            
            % 获取变量名
            tsla_var_name = get_var_name(tsla_data);
            s_var_name = get_var_name(s_data);
            cross_var_name = get_var_name(cross_data);
            
            % 提取数据
            TSLA_AllOrders = tsla_data.(tsla_var_name);
            SSLA_AllOrders = s_data.(s_var_name);
            Cross_AllOrders = cross_data.(cross_var_name);
            
            % 计算HSLA（总展开项）
            HSLA_AllOrders = TSLA_AllOrders + SSLA_AllOrders + Cross_AllOrders;
            
            % 获取数据维度
            [ntime, nlat, nlon, norder] = size(TSLA_AllOrders);
            fprintf('    数据维度: time=%d, lat=%d, lon=%d, order=%d\n', ntime, nlat, nlon, norder);
            
            % 为每阶计算趋势
            fprintf('    计算趋势...\n');
            
            % 初始化趋势数组
            trend_TSLA = zeros(nlat, nlon, norder);
            trend_SSLA = zeros(nlat, nlon, norder);
            trend_HSLA = zeros(nlat, nlon, norder);
            trend_Cross = zeros(nlat, nlon, norder);
            
            % 显著性检验数组
            sig_TSLA = zeros(nlat, nlon, norder);
            sig_SSLA = zeros(nlat, nlon, norder);
            sig_HSLA = zeros(nlat, nlon, norder);
            sig_Cross = zeros(nlat, nlon, norder);
            
            % 进度显示
            fprintf('    进度: ');
            total_points = nlat * nlon * norder;
            processed = 0;
            
            % 对每个网格点和每阶计算趋势
            for order = 1:norder
                for lat_idx = 1:nlat
                    for lon_idx = 1:nlon
                        % TSLA趋势
                        tsla_ts = squeeze(TSLA_AllOrders(:, lat_idx, lon_idx, order));
                        [trend_TSLA(lat_idx, lon_idx, order), sig_TSLA(lat_idx, lon_idx, order)] = ...
                            calculate_ols_trend(time_years(1:ntime), tsla_ts);
                        
                        % SSLA趋势
                        ssla_ts = squeeze(SSLA_AllOrders(:, lat_idx, lon_idx, order));
                        [trend_SSLA(lat_idx, lon_idx, order), sig_SSLA(lat_idx, lon_idx, order)] = ...
                            calculate_ols_trend(time_years(1:ntime), ssla_ts);
                        
                        % HSLA趋势
                        hsla_ts = squeeze(HSLA_AllOrders(:, lat_idx, lon_idx, order));
                        [trend_HSLA(lat_idx, lon_idx, order), sig_HSLA(lat_idx, lon_idx, order)] = ...
                            calculate_ols_trend(time_years(1:ntime), hsla_ts);
                        
                        % Cross趋势
                        cross_ts = squeeze(Cross_AllOrders(:, lat_idx, lon_idx, order));
                        [trend_Cross(lat_idx, lon_idx, order), sig_Cross(lat_idx, lon_idx, order)] = ...
                            calculate_ols_trend(time_years(1:ntime), cross_ts);
                        
                        % 更新进度
                        processed = processed + 1;
                        if mod(processed, round(total_points/10)) == 0
                            fprintf('.');
                        end
                    end
                end
            end
            fprintf(' 完成\n');
            
            % 保存结果
            output_file = fullfile(output_dir, sprintf('%s_%s_Trends.mat', dataset.name, state));
            fprintf('    保存结果到: %s\n', output_file);
            
            save(output_file, 'trend_TSLA', 'trend_SSLA', 'trend_HSLA', 'trend_Cross', ...
                 'sig_TSLA', 'sig_SSLA', 'sig_HSLA', 'sig_Cross', ...
                 'time_years', 'dataset', 'state', 'nlat', 'nlon', 'norder');
            
            % 生成汇总统计
            generate_summary_statistics(dataset.name, state_cn{s_idx}, ...
                trend_TSLA, trend_SSLA, trend_HSLA, trend_Cross, ...
                sig_TSLA, sig_SSLA, sig_HSLA, sig_Cross, output_dir);
            
        catch ME
            fprintf('    错误: %s\n', ME.message);
            continue;
        end
    end
end

fprintf('\n========== 所有趋势计算完成 ==========\n');
fprintf('结果保存在: %s\n', output_dir);

end

function var_name = get_var_name(data_struct)
% 获取MAT文件中的变量名
fields = fieldnames(data_struct);
% 寻找可能的变量名
possible_names = {'TSLA_AllOrders', 'SSLA_AllOrders', 'Cross_AllOrders'};
for i = 1:length(possible_names)
    if isfield(data_struct, possible_names{i})
        var_name = possible_names{i};
        return;
    end
end
% 如果没有找到标准名称，使用第一个非标准字段
var_name = fields{1};
end

function [trend, significance] = calculate_ols_trend(time, data)
% 最小二乘法一元线性回归计算趋势
% trend: 趋势斜率（单位/年）
% significance: 显著性（1表示显著，0表示不显著，p<0.05）

% 移除NaN值
valid_idx = ~isnan(data);
if sum(valid_idx) < 10  % 至少需要10个有效点
    trend = NaN;
    significance = 0;
    return;
end

time_valid = time(valid_idx);
data_valid = data(valid_idx);

% 中心化时间（提高数值稳定性）
time_mean = mean(time_valid);
time_centered = time_valid - time_mean;

% 最小二乘法计算斜率
n = length(time_valid);
Sxy = sum(time_centered .* data_valid);
Sxx = sum(time_centered .^ 2);

if Sxx == 0
    trend = 0;
    significance = 0;
    return;
end

trend = Sxy / Sxx;  % 斜率（趋势）

% 计算截距
intercept = mean(data_valid) - trend * time_mean;

% 计算残差和标准误差
predicted = intercept + trend * time_valid;
residuals = data_valid - predicted;
SSE = sum(residuals .^ 2);
MSE = SSE / (n - 2);  % 均方误差
SE_slope = sqrt(MSE / Sxx);  % 斜率的标准误差

% t检验
if SE_slope == 0
    t_stat = Inf;
else
    t_stat = abs(trend) / SE_slope;
end

% 计算p值（双侧检验）
df = n - 2;  % 自由度
p_value = 2 * (1 - tcdf(t_stat, df));

% 显著性判断（p < 0.05）
significance = p_value < 0.05;

end

function generate_summary_statistics(dataset_name, state_name, ...
trend_TSLA, trend_SSLA, trend_HSLA, trend_Cross, ...
sig_TSLA, sig_SSLA, sig_HSLA, sig_Cross, output_dir)
% 生成汇总统计信息

fprintf('    生成汇总统计...\n');

% 创建汇总结构
summary = struct();
summary.dataset = dataset_name;
summary.state = state_name;
summary.timestamp = datestr(now);

% 计算各阶的全局统计
norder = size(trend_TSLA, 3);
for order = 1:norder
    field_name = sprintf('order%d', order);
    
    % TSLA统计
    tsla_data = trend_TSLA(:,:,order);
    tsla_valid = tsla_data(~isnan(tsla_data));
    summary.(field_name).TSLA.mean = mean(tsla_valid);
    summary.(field_name).TSLA.std = std(tsla_valid);
    summary.(field_name).TSLA.min = min(tsla_valid);
    summary.(field_name).TSLA.max = max(tsla_valid);
    summary.(field_name).TSLA.sig_ratio = mean(sig_TSLA(:,:,order), 'all', 'omitnan');
    
    % SSLA统计
    ssla_data = trend_SSLA(:,:,order);
    ssla_valid = ssla_data(~isnan(ssla_data));
    summary.(field_name).SSLA.mean = mean(ssla_valid);
    summary.(field_name).SSLA.std = std(ssla_valid);
    summary.(field_name).SSLA.min = min(ssla_valid);
    summary.(field_name).SSLA.max = max(ssla_valid);
    summary.(field_name).SSLA.sig_ratio = mean(sig_SSLA(:,:,order), 'all', 'omitnan');
    
    % HSLA统计
    hsla_data = trend_HSLA(:,:,order);
    hsla_valid = hsla_data(~isnan(hsla_data));
    summary.(field_name).HSLA.mean = mean(hsla_valid);
    summary.(field_name).HSLA.std = std(hsla_valid);
    summary.(field_name).HSLA.min = min(hsla_valid);
    summary.(field_name).HSLA.max = max(hsla_valid);
    summary.(field_name).HSLA.sig_ratio = mean(sig_HSLA(:,:,order), 'all', 'omitnan');
    
    % Cross统计
    cross_data = trend_Cross(:,:,order);
    cross_valid = cross_data(~isnan(cross_data));
    summary.(field_name).Cross.mean = mean(cross_valid);
    summary.(field_name).Cross.std = std(cross_valid);
    summary.(field_name).Cross.min = min(cross_valid);
    summary.(field_name).Cross.max = max(cross_valid);
    summary.(field_name).Cross.sig_ratio = mean(sig_Cross(:,:,order), 'all', 'omitnan');
end

% 保存汇总
summary_file = fullfile(output_dir, sprintf('%s_%s_Summary.mat', dataset_name, state_name));
save(summary_file, 'summary');

% 生成文本报告
report_file = fullfile(output_dir, sprintf('%s_%s_Summary.txt', dataset_name, state_name));
fid = fopen(report_file, 'w');
if fid == -1
    fprintf('    警告: 无法创建报告文件 %s\n', report_file);
    return;
end

fprintf(fid, '数据集: %s\n', dataset_name);
fprintf(fid, '状态: %s\n', state_name);
fprintf(fid, '生成时间: %s\n', datestr(now));
fprintf(fid, '\n');
fprintf(fid, '='*60);
fprintf(fid, '\n');
fprintf(fid, '各阶趋势统计（单位：m/年）\n');
fprintf(fid, '='*60);
fprintf(fid, '\n');
fprintf(fid, '阶数  |  TSLA均值    |  SSLA均值    |  HSLA均值    |  Cross均值   |\n');
fprintf(fid, '-----|-------------|-------------|-------------|-------------|\n');

for order = 1:norder
    fprintf(fid, '%d    |  %8.6f  |  %8.6f  |  %8.6f  |  %8.6f  |\n', ...
        order, ...
        summary.(sprintf('order%d', order)).TSLA.mean, ...
        summary.(sprintf('order%d', order)).SSLA.mean, ...
        summary.(sprintf('order%d', order)).HSLA.mean, ...
        summary.(sprintf('order%d', order)).Cross.mean);
end

fprintf(fid, '\n');
fprintf(fid, '='*60);
fprintf(fid, '\n');
fprintf(fid, '显著性比例（p < 0.05）\n');
fprintf(fid, '='*60);
fprintf(fid, '\n');
fprintf(fid, '阶数  |  TSLA显著比例 |  SSLA显著比例 |  HSLA显著比例 |  Cross显著比例 |\n');
fprintf(fid, '-----|--------------|--------------|--------------|--------------|\n');

for order = 1:norder
    fprintf(fid, '%d    |  %8.2f%%    |  %8.2f%%    |  %8.2f%%    |  %8.2f%%    |\n', ...
        order, ...
        summary.(sprintf('order%d', order)).TSLA.sig_ratio * 100, ...
        summary.(sprintf('order%d', order)).SSLA.sig_ratio * 100, ...
        summary.(sprintf('order%d', order)).HSLA.sig_ratio * 100, ...
        summary.(sprintf('order%d', order)).Cross.sig_ratio * 100);
end

fclose(fid);
fprintf('    汇总报告已生成: %s\n', report_file);

end