function Calculate_Hybrid_Trends()
% 用最小二乘法一元线性回归 (OLS Linear Regression) 方法
% 计算三套数据的平均态和标准态的混合项趋势文件
% 整合了：数据维度检测、单位智能统一、极速向量化计算(矩阵左除)

fprintf('========== 计算混合项趋势（极速增强版） ==========\n\n');

% 基础目录
base_dir = 'D:\work';
output_dir = fullfile(base_dir, 'Task_Convergence', 'Trend_Results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 三套数据目录
data_dirs = {
    struct('name', 'EN4', 'path', fullfile(base_dir, 'EN4_TSLA_Terms'));
    struct('name', 'IAP', 'path', fullfile(base_dir, 'IAP_mat_data'));
    struct('name', 'Ishii', 'path', fullfile(base_dir, 'Ishii_mat_data'));
};

% 两种状态
states = {'Average', 'StdRef'};
state_cn = {'平均态', '标准态'};

% 主循环：计算每套数据、每种状态、每阶的趋势
for d_idx = 1:length(data_dirs)
    dataset = data_dirs{d_idx};
    fprintf('\n处理数据集: %s\n', dataset.name);
    
    for s_idx = 1:length(states)
        state = states{s_idx};
        fprintf('  状态: %s (%s)\n', state_cn{s_idx}, state);
        
        % 构建文件名（尝试多种可能的格式）
        cross_file = [];
        
        % 尝试第一种格式: {dataset}_Cross_Terms_1to8_{state}.mat
        cross_file1 = fullfile(dataset.path, sprintf('%s_Cross_Terms_1to8_%s.mat', dataset.name, state));
        
        % 尝试第二种格式: {dataset}_Cross_Terms_{state}.mat
        cross_file2 = fullfile(dataset.path, sprintf('%s_Cross_Terms_%s.mat', dataset.name, state));
        
        % 尝试第三种格式: Cross_Terms_1to8_{state}.mat
        cross_file3 = fullfile(dataset.path, sprintf('Cross_Terms_1to8_%s.mat', state));
        
        % 尝试第四种格式: Cross_Terms_{state}.mat
        cross_file4 = fullfile(dataset.path, sprintf('Cross_Terms_%s.mat', state));
        
        % 选择存在的文件
        if exist(cross_file1, 'file'), cross_file = cross_file1; end
        if exist(cross_file2, 'file'), cross_file = cross_file2; end
        if exist(cross_file3, 'file'), cross_file = cross_file3; end
        if exist(cross_file4, 'file'), cross_file = cross_file4; end
        
        % 检查输入文件是否存在
        if isempty(cross_file)
            fprintf('    警告: 缺少 Cross 泰勒展开文件\n');
            continue;
        end
        
        % 显示找到的文件
        fprintf('    找到文件:\n');
        fprintf('    - Cross: %s\n', cross_file);
        
        % 检查输出文件是否已存在
        output_file = fullfile(output_dir, sprintf('%s_%s_Hybrid_Trends.mat', dataset.name, state));
        
        if exist(output_file, 'file')
            fprintf('    跳过: 混合项趋势文件已存在 %s\n', output_file);
            continue;
        end
        
        fprintf('    开始处理...\n')
        
        % 加载数据
        fprintf('    加载数据...\n');
        try
            % 加载Cross数据
            cross_data = load(cross_file);
            
            % 检查数据结构
            Check_Data_Structure(cross_data, cross_file);
            
            % 获取时间轴
            if isfield(cross_data, 'Time_Axis')
                time_vec = cross_data.Time_Axis(:);
            elseif isfield(cross_data, 'time_axis')
                time_vec = cross_data.time_axis(:);
            elseif isfield(cross_data, 'time_vec')
                time_vec = cross_data.time_vec(:);
            else
                % 如果没有时间轴，使用默认时间（2005-2024年，每月）
                fprintf('    警告: 未找到时间轴，使用默认时间向量\n');
                years = 2005:2024;
                months = 1:12;
                time_vec = zeros(length(years) * length(months), 1);
                for i = 1:length(years)
                    for j = 1:length(months)
                        idx = (i-1)*12 + j;
                        time_vec(idx) = years(i) + (months(j)-1)/12;
                    end
                end
            end
            
            % 获取坐标
            if isfield(cross_data, 'Lon')
                Lon = cross_data.Lon;
            elseif isfield(cross_data, 'lon')
                Lon = cross_data.lon;
            else
                Lon = [];
                fprintf('    警告: 未找到经度坐标\n');
            end
            
            if isfield(cross_data, 'Lat')
                Lat = cross_data.Lat;
            elseif isfield(cross_data, 'lat')
                Lat = cross_data.lat;
            else
                Lat = [];
                fprintf('    警告: 未找到纬度坐标\n');
            end
            
            % 获取变量名
            cross_var_name = get_var_name(cross_data, 'Cross_AllOrders');
            Cross_AllOrders = cross_data.(cross_var_name);
            
            % 检测数据维度顺序
            [dim1, dim2, dim3, dim4] = size(Cross_AllOrders);
            fprintf('    原始数据维度: %d x %d x %d x %d\n', dim1, dim2, dim3, dim4);
            
            % 确定时间维度位置
            time_dim = 0;
            if dim1 == length(time_vec)
                time_dim = 1;
                fprintf('    时间维度在第1维\n');
            elseif dim3 == length(time_vec)
                time_dim = 3;
                fprintf('    时间维度在第3维\n');
            else
                fprintf('    警告: 时间维度长度不匹配，使用默认处理\n');
                time_dim = 1;
            end
            
            % --- 单位统一处理 ---
            % 如果数据极值小于10，认定为"米"，强制转为"毫米"
            max_cross = max(abs(Cross_AllOrders(:)));
            if ~isnan(max_cross) && max_cross < 10
                fprintf('    [自动修正] 发现数据单位疑似为"米"，已自动乘以1000转换为"毫米"\n');
                Cross_AllOrders = Cross_AllOrders * 1000;
            else
                fprintf('    [单位确认] 数据单位已确认，无需转换\n');
            end
            
            % 获取并规范数据维度 (最终全统一为: [Time, Lat, Lon, Order])
            if time_dim == 1
                [ntime, nlat, nlon, norder] = size(Cross_AllOrders);
            else
                [nlon, nlat, ntime, norder] = size(Cross_AllOrders);
                % 重新排列维度为 [Time, Lat, Lon, Order]
                Cross_AllOrders = permute(Cross_AllOrders, [3, 2, 1, 4]);
                % 调整坐标顺序
                if ~isempty(Lon) && size(Lon,1) > 1
                    Lon = Lon';
                end
                if ~isempty(Lat) && size(Lat,2) > 1
                    Lat = Lat';
                end
            end
            
            fprintf('    处理后数据维度: time=%d, lat=%d, lon=%d, order=%d\n', ntime, nlat, nlon, norder);
            
            % =========================================================
            % 🚀 极速向量化计算趋势与显著性 (无缝替换原来的三重循环)
            % =========================================================
            fprintf('    极速向量化计算混合项趋势与显著性...\n');
            
            % 初始化趋势数组
            trend_Hybrid = zeros(nlat, nlon, norder);
            
            % 显著性检验数组
            sig_Hybrid = zeros(nlat, nlon, norder);
            
            % 构建公共回归矩阵 X
            time_valid = time_vec(1:ntime);
            
            % 中心化时间轴，提高数值稳定性
            time_centered = time_valid - mean(time_valid);
            X = [time_centered, ones(ntime, 1)];
            
            % 计算Sxx和自由度
            Sxx = sum(time_centered.^2);
            df = ntime - 2;
            
            % 检查X矩阵的秩，避免秩亏
            rank_X = rank(X);
            if rank_X < size(X, 2)
                fprintf('    [警告] 时间轴数据存在问题，X矩阵秩亏 (rank=%d)，将使用简化回归\n', rank_X);
                % 使用简化的回归方法
                use_simple_regression = true;
            else
                use_simple_regression = false;
            end
            
            % 对每阶进行向量化计算
            for order = 1:norder
                % 将三维切片压扁为二维: (Time, Lat*Lon)
                Y_Hybrid = reshape(Cross_AllOrders(:,:,:,order), ntime, nlat*nlon);
                
                % 找到非NaN的有效海洋点（只校验第一个时间步）
                valid_idx = ~isnan(Y_Hybrid(1,:));
                
                % --- 混合项 ---                
                [t_Hybrid, s_Hybrid] = batch_ols(X, Y_Hybrid(:, valid_idx), Sxx, df);
                trend_1d = nan(1, nlat*nlon); sig_1d = zeros(1, nlat*nlon);
                trend_1d(valid_idx) = t_Hybrid; sig_1d(valid_idx) = s_Hybrid;
                trend_Hybrid(:,:,order) = reshape(trend_1d, nlat, nlon);
                sig_Hybrid(:,:,order) = reshape(sig_1d, nlat, nlon);
                
                fprintf('      第 %d 阶计算完成\n', order);
            end
            fprintf('    极速计算全部完成！\n');
            % =========================================================
            
            % 保存结果
            output_file = fullfile(output_dir, sprintf('%s_%s_Hybrid_Trends.mat', dataset.name, state));
            fprintf('    保存结果到: %s\n', output_file);
            
            save(output_file, 'trend_Hybrid', 'sig_Hybrid', ...
                 'time_vec', 'Lon', 'Lat', 'dataset', 'state', 'nlat', 'nlon', 'norder');
            
            % 生成汇总统计
            generate_hybrid_summary_statistics(dataset.name, state_cn{s_idx}, ...
                trend_Hybrid, sig_Hybrid, output_dir);
            
        catch ME
            fprintf('    错误: %s\n', ME.message);
            continue;
        end
    end
end

fprintf('\n========== 所有混合项趋势计算完成 ==========\n');
fprintf('结果保存在: %s\n', output_dir);

end

% =========================================================================
% 极速批量 OLS 趋势与显著性计算函数
% =========================================================================
function [trends, sigs] = batch_ols(X, Y_valid, Sxx, df)
    if isempty(Y_valid)
        trends = []; sigs = []; return;
    end
    
    % 检查X矩阵是否秩亏
    rank_X = rank(X);
    if rank_X < size(X, 2)
        % 秩亏情况：使用简化的回归方法
        fprintf('    [处理] X矩阵秩亏，使用简化回归方法\n');
        
        % 只使用中心化时间轴进行简单线性回归
        time_centered = X(:, 1);
        n_points = size(time_centered, 1);
        
        % 对每个网格点单独计算
        n_points = size(Y_valid, 2);
        trends = zeros(1, n_points);
        sigs = zeros(1, n_points);
        
        for i = 1:n_points
            y = Y_valid(:, i);
            valid_y = y(~isnan(y));
            valid_t = time_centered(~isnan(y));
            
            if length(valid_y) >= 10
                % 简单线性回归
                Sxy = sum(valid_t .* valid_y);
                if Sxx > 1e-10
                    trend = Sxy / Sxx;
                    trends(i) = trend;
                    
                    % 计算显著性
                    y_mean = mean(valid_y);
                    y_pred = trend * valid_t + y_mean;
                    res = valid_y - y_pred;
                    mse = sum(res.^2) / (length(valid_y) - 2);
                    se = sqrt(mse / Sxx);
                    if se > 0
                        t_stat = abs(trend) / se;
                        p_value = 2 * (1 - tcdf(t_stat, length(valid_y) - 2));
                        sigs(i) = p_value < 0.05;
                    end
                end
            end
        end
    else
        % 正常情况：使用矩阵左除
        B = X \ Y_valid; 
        trends = B(1, :);
        
        % 批量显著性 t 检验
        Y_pred = X * B;
        Res = Y_valid - Y_pred;
        SSE = sum(Res.^2, 1);
        MSE = SSE / df;
        SE_slope = sqrt(MSE / Sxx);
        
        t_stat = abs(trends) ./ SE_slope;
        
        % tcdf 函数求 p 值（需 Statistics and Machine Learning Toolbox）
        p_values = 2 * (1 - tcdf(t_stat, df)); 
        sigs = p_values < 0.05;
    end
end

% =========================================================================
% 下面是各类辅助检查函数（获取变量名、维度体检、报告生成）
% =========================================================================

function var_name = get_var_name(data_struct, default_name)
    fields = fieldnames(data_struct);
    
    % 1. 精确匹配默认名称
    if isfield(data_struct, default_name)
        var_name = default_name;
        return;
    end
    
    % 2. 模糊匹配常用的包含词
    for i = 1:length(fields)
        if contains(fields{i}, 'AllOrders') || contains(fields{i}, 'Cross')
            var_name = fields{i};
            return;
        end
    end
    
    % 3. 排除已知坐标变量名
    exclude_vars = {'lon', 'lat', 'time', 'Lon', 'Lat', 'Time_Axis', 'time_vec', 'time_axis'};
    for i = 1:length(fields)
        if ~ismember(fields{i}, exclude_vars)
            var_name = fields{i};
            return;
        end
    end
    
    var_name = fields{1}; % 兜底
end

function Check_Data_Structure(data, filepath) 
     fprintf('\n正在检查数据体检报告: %s\n', filepath); 
     
     vars = fieldnames(data); 
     var_idx = find(contains(vars, 'AllOrders')); 
     if isempty(var_idx) 
         error('文件中未找到包含 "AllOrders" 的变量！请检查数据文件。'); 
     end 
     var_name = vars{var_idx(1)}; 
     ts_data = data.(var_name); 
     
     fprintf('========== 数据体检报告: %s ==========\n', var_name); 
     
     dims = size(ts_data); 
     fprintf('1. 矩阵维度: ['); fprintf('%d ', dims); fprintf(']\n'); 
     
     if length(dims) == 4 
         slice = ts_data(:, :, 1, 1); 
     else 
         slice = ts_data(:, :, 1); 
     end 
     
     num_nans = sum(isnan(slice(:))); 
     num_zeros = sum(slice(:) == 0); 
     total_pts = numel(slice); 
     
     fprintf('2. 陆地/缺测值检查:\n'); 
     fprintf('   -> NaN 数量: %d (占比 %.1f%%)\n', num_nans, num_nans/total_pts*100); 
     fprintf('   -> 绝对0值数量: %d (占比 %.1f%%)\n', num_zeros, num_zeros/total_pts*100); 
     if num_zeros > total_pts * 0.1 
         fprintf('   [!!! 警告 !!!] 发现大量 0 值，陆地可能被错误填成了 0！\n'); 
     end 
     
     valid_data = ts_data(~isnan(ts_data)); 
     max_val = max(valid_data(:)); 
     min_val = min(valid_data(:)); 
     fprintf('3. 数值量级检查:\n'); 
     fprintf('   -> 整体最大值: %.4f\n', max_val); 
     fprintf('   -> 整体最小值: %.4f\n', min_val); 
     
     fprintf('4. 数据单位检查 (Unit Check):\n'); 
     abs_max = max(abs(max_val), abs(min_val)); 
     if abs_max > 10 && abs_max < 5000 
         fprintf('   [推断结果] 数据极值为 %.2f，量级在百位级别。\n', abs_max); 
         fprintf('   [推断单位] => 毫米 (mm) (趋势单位将是 mm/yr)\n'); 
     elseif abs_max <= 10 
         fprintf('   [推断结果] 数据极值为 %.2f，量级非常小 (<10)。\n', abs_max); 
         fprintf('   [!!! 严重警告 !!!] => 该数据极大概率是以 米 (m) 为单位！\n'); 
     else 
         fprintf('   [推断结果] 数据极值异常大 (%.2f)，请核查数据是否正确！\n', abs_max); 
     end 
     fprintf('==================================================\n\n'); 
end

function generate_hybrid_summary_statistics(dataset_name, state_name, ...
trend_Hybrid, sig_Hybrid, output_dir)

fprintf('    生成混合项汇总统计...\n');
summary = struct();
summary.dataset = dataset_name;
summary.state = state_name;
summary.timestamp = datestr(now);

norder = size(trend_Hybrid, 3);
for order = 1:norder
    field_name = sprintf('order%d', order);
    
    hybrid_data = trend_Hybrid(:,:,order);
    hybrid_valid = hybrid_data(~isnan(hybrid_data));
    summary.(field_name).Hybrid.mean = mean(hybrid_valid);
    summary.(field_name).Hybrid.std = std(hybrid_valid);
    summary.(field_name).Hybrid.min = min(hybrid_valid);
    summary.(field_name).Hybrid.max = max(hybrid_valid);
    summary.(field_name).Hybrid.sig_ratio = mean(sig_Hybrid(:,:,order), 'all', 'omitnan');
end

summary_file = fullfile(output_dir, sprintf('%s_%s_Hybrid_Summary.mat', dataset_name, state_name));
save(summary_file, 'summary');

report_file = fullfile(output_dir, sprintf('%s_%s_Hybrid_Summary.txt', dataset_name, state_name));
fid = fopen(report_file, 'w');
if fid == -1, return; end

fprintf(fid, '数据集: %s\n状态: %s\n生成时间: %s\n\n', dataset_name, state_name, datestr(now));
fprintf(fid, '============================================================\n');
fprintf(fid, '各阶混合项趋势统计（单位：mm/年）\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '阶数  |  混合项均值    |  混合项标准差  |  混合项最小值  |  混合项最大值  |\n');
fprintf(fid, '-----|-------------|-------------|-------------|-------------|\n');

for order = 1:norder
    hybrid_m = summary.(sprintf('order%d', order)).Hybrid.mean;
    hybrid_std = summary.(sprintf('order%d', order)).Hybrid.std;
    hybrid_min = summary.(sprintf('order%d', order)).Hybrid.min;
    hybrid_max = summary.(sprintf('order%d', order)).Hybrid.max;
    fprintf(fid, '%d    |  %8.6f  |  %8.6f  |  %8.6f  |  %8.6f  |\n', order, hybrid_m, hybrid_std, hybrid_min, hybrid_max);
end

fprintf(fid, '\n============================================================\n');
fprintf(fid, '显著性比例（p < 0.05）\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '阶数  |  混合项显著比例 |\n');
fprintf(fid, '-----|--------------|\n');

for order = 1:norder
    hybrid_s = summary.(sprintf('order%d', order)).Hybrid.sig_ratio * 100;
    fprintf(fid, '%d    |  %8.2f%%    |\n', order, hybrid_s);
end
fclose(fid);
fprintf('    混合项汇总报告已生成: %s\n', report_file);
end