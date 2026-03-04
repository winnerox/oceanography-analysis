function Calculate_Amplitude_Enhanced()
% 用标准差法 (Standard Deviation / STD) 方法
% 计算三套数据的平均态和标准态的振幅文件
% 整合了：数据维度检测、单位智能统一、极速向量化计算

fprintf('========== 计算泰勒展开项振幅（极速增强版） ==========\n\n');

% 基础目录
base_dir = 'D:\work';
output_dir = fullfile(base_dir, 'Task_Convergence', 'Amplitude_Results');
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

% 主循环：计算每套数据、每种状态、每阶的振幅
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
        
        % 检查输入文件是否存在
        if ~exist(tsla_file, 'file')
            fprintf('    警告: 文件不存在 %s\n', tsla_file);
            continue;
        end
        
        % 检查输出文件是否已存在
        output_file = fullfile(output_dir, sprintf('%s_%s_Amplitude.mat', dataset.name, state));
        gemini_output_file = fullfile(dataset.path, sprintf('%s_Amplitude_1to8_%s.mat', dataset.name, state));
        
        if exist(output_file, 'file') && exist(gemini_output_file, 'file')
            fprintf('    跳过: 振幅文件已存在 %s\n', output_file);
            continue;
        end
        
        fprintf('    开始处理...\n')
        
        % 加载数据
        fprintf('    加载数据...\n');
        try
            % 加载TSLA数据和坐标信息
            tsla_data = load(tsla_file);
            
            % 检查数据结构
            Check_Data_Structure(tsla_file);
            
            % 获取时间轴
            if isfield(tsla_data, 'Time_Axis')
                time_vec = tsla_data.Time_Axis(:);
            elseif isfield(tsla_data, 'time_axis')
                time_vec = tsla_data.time_axis(:);
            elseif isfield(tsla_data, 'time_vec')
                time_vec = tsla_data.time_vec(:);
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
            if isfield(tsla_data, 'Lon')
                Lon = tsla_data.Lon;
            elseif isfield(tsla_data, 'lon')
                Lon = tsla_data.lon;
            else
                Lon = [];
                fprintf('    警告: 未找到经度坐标\n');
            end
            
            if isfield(tsla_data, 'Lat')
                Lat = tsla_data.Lat;
            elseif isfield(tsla_data, 'lat')
                Lat = tsla_data.lat;
            else
                Lat = [];
                fprintf('    警告: 未找到纬度坐标\n');
            end
            
            % 获取变量名
            tsla_var_name = get_var_name(tsla_data, 'TSLA_AllOrders');
            TSLA_AllOrders = tsla_data.(tsla_var_name);
            
            % 检测数据维度顺序
            [dim1, dim2, dim3, dim4] = size(TSLA_AllOrders);
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
            
            % 加载SSLA数据
            SSLA_AllOrders = [];
            if exist(s_file, 'file')
                s_data = load(s_file);
                s_var_name = get_var_name(s_data, 'SSLA_AllOrders');
                SSLA_AllOrders = s_data.(s_var_name);
            end
            
            % 加载Cross数据
            Cross_AllOrders = [];
            if exist(cross_file, 'file')
                cross_data = load(cross_file);
                cross_var_name = get_var_name(cross_data, 'Cross_AllOrders');
                Cross_AllOrders = cross_data.(cross_var_name);
            end
            
            % 计算HSLA（总展开项）
            if ~isempty(SSLA_AllOrders) && ~isempty(Cross_AllOrders)
                HSLA_AllOrders = TSLA_AllOrders + SSLA_AllOrders + Cross_AllOrders;
            else
                HSLA_AllOrders = [];
            end
            
            % --- 单位统一处理 ---
            % 如果数据极值小于10，认定为"米"，强制转为"毫米"
            max_tsla = max(abs(TSLA_AllOrders(:)));
            if ~isnan(max_tsla) && max_tsla < 10
                fprintf('    [自动修正] 发现数据单位疑似为"米"，已自动乘以1000转换为"毫米"\n');
                TSLA_AllOrders = TSLA_AllOrders * 1000;
                if ~isempty(SSLA_AllOrders), SSLA_AllOrders = SSLA_AllOrders * 1000; end
                if ~isempty(Cross_AllOrders), Cross_AllOrders = Cross_AllOrders * 1000; end
                if ~isempty(HSLA_AllOrders), HSLA_AllOrders = HSLA_AllOrders * 1000; end
            else
                fprintf('    [单位确认] 数据单位已确认，无需转换\n');
            end
            
            % 获取并规范数据维度 (最终全统一为: [Time, Lat, Lon, Order])
            if time_dim == 1
                [ntime, nlat, nlon, norder] = size(TSLA_AllOrders);
            else
                [nlon, nlat, ntime, norder] = size(TSLA_AllOrders);
                % 重新排列维度为 [Time, Lat, Lon, Order]
                TSLA_AllOrders = permute(TSLA_AllOrders, [3, 2, 1, 4]);
                if ~isempty(SSLA_AllOrders)
                    SSLA_AllOrders = permute(SSLA_AllOrders, [3, 2, 1, 4]);
                end
                if ~isempty(Cross_AllOrders)
                    Cross_AllOrders = permute(Cross_AllOrders, [3, 2, 1, 4]);
                end
                if ~isempty(HSLA_AllOrders)
                    HSLA_AllOrders = permute(HSLA_AllOrders, [3, 2, 1, 4]);
                end
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
            % 🚀 极速向量化计算振幅 (去趋势残差标准差法)
            % =========================================================
            fprintf('    极速向量化计算振幅...\n');
            
            % 初始化振幅数组
            amp_TSLA = zeros(nlat, nlon, norder);
            amp_SSLA = zeros(nlat, nlon, norder);
            amp_HSLA = zeros(nlat, nlon, norder);
            amp_Cross = zeros(nlat, nlon, norder);
            
            % 准备时间向量用于去趋势
            time_valid = time_vec(1:ntime);
            time_centered = time_valid - mean(time_valid);
            X_reg = [time_centered, ones(ntime, 1)];
            
            % 对每阶进行向量化计算
            for order = 1:norder
                % --- TSLA ---
                % 将三维切片压扁为二维: (Time, Lat*Lon)
                Y_TSLA = reshape(TSLA_AllOrders(:,:,:,order), ntime, nlat*nlon);
                
                % 找到非NaN的有效海洋点
                valid_idx = ~isnan(Y_TSLA(1,:));
                Y_TSLA_valid = Y_TSLA(:, valid_idx);
                
                % 去趋势计算
                if size(Y_TSLA_valid, 2) > 0
                    % 最小二乘拟合
                    B = X_reg \ Y_TSLA_valid;
                    % 计算残差
                    Res_TSLA = Y_TSLA_valid - X_reg * B;
                    % 计算残差的标准差
                    amp_TSLA_flat = nanstd(Res_TSLA, 0, 1); % 0 表示除以 n-1
                    
                    % 填充结果
                    amp_TSLA_flat_full = nan(1, nlat*nlon);
                    amp_TSLA_flat_full(valid_idx) = amp_TSLA_flat;
                    amp_TSLA(:,:,order) = reshape(amp_TSLA_flat_full, nlat, nlon);
                end
                
                % --- SSLA ---
                if ~isempty(SSLA_AllOrders)
                    Y_SSLA = reshape(SSLA_AllOrders(:,:,:,order), ntime, nlat*nlon);
                    Y_SSLA_valid = Y_SSLA(:, valid_idx);
                    
                    if size(Y_SSLA_valid, 2) > 0
                        B = X_reg \ Y_SSLA_valid;
                        Res_SSLA = Y_SSLA_valid - X_reg * B;
                        amp_SSLA_flat = nanstd(Res_SSLA, 0, 1);
                        
                        amp_SSLA_flat_full = nan(1, nlat*nlon);
                        amp_SSLA_flat_full(valid_idx) = amp_SSLA_flat;
                        amp_SSLA(:,:,order) = reshape(amp_SSLA_flat_full, nlat, nlon);
                    end
                end
                
                % --- HSLA ---
                if ~isempty(HSLA_AllOrders)
                    Y_HSLA = reshape(HSLA_AllOrders(:,:,:,order), ntime, nlat*nlon);
                    Y_HSLA_valid = Y_HSLA(:, valid_idx);
                    
                    if size(Y_HSLA_valid, 2) > 0
                        B = X_reg \ Y_HSLA_valid;
                        Res_HSLA = Y_HSLA_valid - X_reg * B;
                        amp_HSLA_flat = nanstd(Res_HSLA, 0, 1);
                        
                        amp_HSLA_flat_full = nan(1, nlat*nlon);
                        amp_HSLA_flat_full(valid_idx) = amp_HSLA_flat;
                        amp_HSLA(:,:,order) = reshape(amp_HSLA_flat_full, nlat, nlon);
                    end
                end
                
                % --- Cross ---
                if ~isempty(Cross_AllOrders)
                    Y_Cross = reshape(Cross_AllOrders(:,:,:,order), ntime, nlat*nlon);
                    Y_Cross_valid = Y_Cross(:, valid_idx);
                    
                    if size(Y_Cross_valid, 2) > 0
                        B = X_reg \ Y_Cross_valid;
                        Res_Cross = Y_Cross_valid - X_reg * B;
                        amp_Cross_flat = nanstd(Res_Cross, 0, 1);
                        
                        amp_Cross_flat_full = nan(1, nlat*nlon);
                        amp_Cross_flat_full(valid_idx) = amp_Cross_flat;
                        amp_Cross(:,:,order) = reshape(amp_Cross_flat_full, nlat, nlon);
                    end
                end
                
                fprintf('      第 %d 阶计算完成\n', order);
            end
            fprintf('    极速计算全部完成！\n');
            % =========================================================
            
            % 保存结果
            output_file = fullfile(output_dir, sprintf('%s_%s_Amplitude.mat', dataset.name, state));
            fprintf('    保存结果到: %s\n', output_file);
            
            save(output_file, 'amp_TSLA', 'amp_SSLA', 'amp_HSLA', 'amp_Cross', ...
                 'time_vec', 'Lon', 'Lat', 'dataset', 'state', 'nlat', 'nlon', 'norder');
            
            % 同时保存到数据集目录（兼容Gemini的保存位置）
            gemini_output_file = fullfile(dataset.path, sprintf('%s_Amplitude_1to8_%s.mat', dataset.name, state));
            fprintf('    同时保存到: %s\n', gemini_output_file);
            save(gemini_output_file, 'amp_TSLA', 'amp_SSLA', 'amp_Cross', 'Lon', 'Lat', '-v7.3');
            
            % 生成汇总统计
            generate_summary_statistics(dataset.name, state_cn{s_idx}, ...
                amp_TSLA, amp_SSLA, amp_HSLA, amp_Cross, output_dir);
            
        catch ME
            fprintf('    错误: %s\n', ME.message);
            continue;
        end
    end
end

fprintf('\n========== 所有振幅计算完成 ==========\n');
fprintf('结果保存在: %s\n', output_dir);

end

% =========================================================================
% 下面是各类辅助检查函数（获取变量名、维度体检、报告生成）
% =========================================================================

function var_name = get_var_name(data_struct, default_name)
    fields = fieldnames(data_struct);
    possible_names = {default_name};
    for i = 1:length(possible_names)
        if isfield(data_struct, possible_names{i})
            var_name = possible_names{i};
            return;
        end
    end
    var_name = fields{1};
end

function Check_Data_Structure(filepath) 
     fprintf('\n正在读取文件: %s\n', filepath); 
     data = load(filepath); 
     
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
         fprintf('   [推断单位] => 毫米 (mm)\n'); 
     elseif abs_max <= 10 
         fprintf('   [推断结果] 数据极值为 %.2f，量级非常小 (<10)。\n', abs_max); 
         fprintf('   [!!! 严重警告 !!!] => 该数据极大概率是以 米 (m) 为单位！\n'); 
     else 
         fprintf('   [推断结果] 数据极值异常大 (%.2f)，请核查数据是否正确！\n', abs_max); 
     end 
     fprintf('==================================================\n\n'); 
end

function generate_summary_statistics(dataset_name, state_name, ...
amp_TSLA, amp_SSLA, amp_HSLA, amp_Cross, output_dir)

fprintf('    生成汇总统计...\n');
summary = struct();
summary.dataset = dataset_name;
summary.state = state_name;
summary.timestamp = datestr(now);

norder = size(amp_TSLA, 3);
for order = 1:norder
    field_name = sprintf('order%d', order);
    
    tsla_data = amp_TSLA(:,:,order);
    tsla_valid = tsla_data(~isnan(tsla_data));
    summary.(field_name).TSLA.mean = mean(tsla_valid);
    summary.(field_name).TSLA.std = std(tsla_valid);
    summary.(field_name).TSLA.min = min(tsla_valid);
    summary.(field_name).TSLA.max = max(tsla_valid);
    
    if ~isempty(amp_SSLA)
        ssla_data = amp_SSLA(:,:,order);
        ssla_valid = ssla_data(~isnan(ssla_data));
        summary.(field_name).SSLA.mean = mean(ssla_valid);
        summary.(field_name).SSLA.std = std(ssla_valid);
        summary.(field_name).SSLA.min = min(ssla_valid);
        summary.(field_name).SSLA.max = max(ssla_valid);
    end
    
    if ~isempty(amp_HSLA)
        hsla_data = amp_HSLA(:,:,order);
        hsla_valid = hsla_data(~isnan(hsla_data));
        summary.(field_name).HSLA.mean = mean(hsla_valid);
        summary.(field_name).HSLA.std = std(hsla_valid);
        summary.(field_name).HSLA.min = min(hsla_valid);
        summary.(field_name).HSLA.max = max(hsla_valid);
    end
    
    if ~isempty(amp_Cross)
        cross_data = amp_Cross(:,:,order);
        cross_valid = cross_data(~isnan(cross_data));
        summary.(field_name).Cross.mean = mean(cross_valid);
        summary.(field_name).Cross.std = std(cross_valid);
        summary.(field_name).Cross.min = min(cross_valid);
        summary.(field_name).Cross.max = max(cross_valid);
    end
end

summary_file = fullfile(output_dir, sprintf('%s_%s_Summary.mat', dataset_name, state_name));
save(summary_file, 'summary');

report_file = fullfile(output_dir, sprintf('%s_%s_Summary.txt', dataset_name, state_name));
fid = fopen(report_file, 'w');
if fid == -1, return; end

fprintf(fid, '数据集: %s\n状态: %s\n生成时间: %s\n\n', dataset_name, state_name, datestr(now));
fprintf(fid, '============================================================\n');
fprintf(fid, '各阶振幅统计（单位：mm）\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '阶数  |  TSLA均值    |  SSLA均值    |  HSLA均值    |  Cross均值   |\n');
fprintf(fid, '-----|-------------|-------------|-------------|-------------|\n');

for order = 1:norder
    tsla_m = summary.(sprintf('order%d', order)).TSLA.mean;
    ssla_m = 0; hsla_m = 0; cross_m = 0;
    if isfield(summary.(sprintf('order%d', order)), 'SSLA'), ssla_m = summary.(sprintf('order%d', order)).SSLA.mean; end
    if isfield(summary.(sprintf('order%d', order)), 'HSLA'), hsla_m = summary.(sprintf('order%d', order)).HSLA.mean; end
    if isfield(summary.(sprintf('order%d', order)), 'Cross'), cross_m = summary.(sprintf('order%d', order)).Cross.mean; end
    fprintf(fid, '%d    |  %8.6f  |  %8.6f  |  %8.6f  |  %8.6f  |\n', order, tsla_m, ssla_m, hsla_m, cross_m);
end

fprintf(fid, '\n============================================================\n');
fprintf(fid, '振幅范围\n');
fprintf(fid, '============================================================\n');
fprintf(fid, '阶数  |  TSLA范围    |  SSLA范围    |  HSLA范围    |  Cross范围   |\n');
fprintf(fid, '-----|-------------|-------------|-------------|-------------|\n');

for order = 1:norder
    tsla_min = summary.(sprintf('order%d', order)).TSLA.min;
    tsla_max = summary.(sprintf('order%d', order)).TSLA.max;
    ssla_min = 0; ssla_max = 0; hsla_min = 0; hsla_max = 0; cross_min = 0; cross_max = 0;
    if isfield(summary.(sprintf('order%d', order)), 'SSLA')
        ssla_min = summary.(sprintf('order%d', order)).SSLA.min;
        ssla_max = summary.(sprintf('order%d', order)).SSLA.max;
    end
    if isfield(summary.(sprintf('order%d', order)), 'HSLA')
        hsla_min = summary.(sprintf('order%d', order)).HSLA.min;
        hsla_max = summary.(sprintf('order%d', order)).HSLA.max;
    end
    if isfield(summary.(sprintf('order%d', order)), 'Cross')
        cross_min = summary.(sprintf('order%d', order)).Cross.min;
        cross_max = summary.(sprintf('order%d', order)).Cross.max;
    end
    fprintf(fid, '%d    |  %.4f-%.4f |  %.4f-%.4f |  %.4f-%.4f |  %.4f-%.4f |\n', ...
        order, tsla_min, tsla_max, ssla_min, ssla_max, hsla_min, hsla_max, cross_min, cross_max);
end
fclose(fid);
fprintf('    汇总报告已生成: %s\n', report_file);
end