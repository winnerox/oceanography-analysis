function Calculate_CrossDetail_Trends()
% =========================================================================
% 功能：极速计算 3 套产品 (EN4, IAP, Ishii) 
%       × 2 种状态 (Average, StdRef) 
%       × 28 个分解混合项 (T1S1, T2S1, ..., T1S7) 的线性空间趋势
%
% 依赖：D:\work\[产品名]_TSLA_Terms\[产品名]_CrossDetail_[态].mat
% 输出：D:\work\Task_Convergence\Trend_Results\[产品名]_[态]_CrossDetail_Trends.mat
% =========================================================================

clc;
fprintf('========== 开始批量计算 28项分解混合项趋势 (CrossDetail) ==========\n\n');

base_dir = 'D:\work';
output_dir = fullfile(base_dir, 'Task_Convergence', 'Trend_Results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 预定义28项变量名列表 (2-8阶的所有S>=1混合项)
var_list = {};
for n = 2:8
    for k = 1:n-1
        var_list{end+1} = sprintf('Cross_T%dS%d', n-k, k);
    end
end
num_vars = length(var_list);

% 三套数据配置
datasets = {
    struct('name', 'EN4',   'folder', fullfile(base_dir, 'EN4_TSLA_Terms'));
    struct('name', 'IAP',   'folder', fullfile(base_dir, 'IAP_mat_data'));
    struct('name', 'Ishii', 'folder', fullfile(base_dir, 'Ishii_mat_data'));
};

% 状态配置
states_in  = {'Avg', 'Std'};
states_out = {'Average', 'StdRef'};
states_cn  = {'平均态', '标准态'};

% 开始大循环
for d = 1:length(datasets)
    ds = datasets{d};
    fprintf('处理数据集: %s\n', ds.name);
    
    for s = 1:length(states_in)
        state_in = states_in{s};
        state_out = states_out{s};
        fprintf('  状态: %s (%s)\n', states_cn{s}, state_out);
        
        % 尝试定位文件
        file_path = fullfile(ds.folder, sprintf('%s_CrossDetail_%s.mat', ds.name, state_in));
        if ~exist(file_path, 'file')
            % 尝试替换可能的命名规则
            file_path2 = fullfile(ds.folder, sprintf('CrossDetail_%s.mat', state_in));
            if exist(file_path2, 'file')
                file_path = file_path2;
            else
                fprintf('    [跳过] 找不到输入文件: %s\n', file_path);
                continue;
            end
        end
        
        % 输出文件
        out_file = fullfile(output_dir, sprintf('%s_%s_CrossDetail_Trends.mat', ds.name, state_out));
        if exist(out_file, 'file')
            fprintf('    [跳过] 趋势文件已存在: %s\n', out_file);
            continue;
        end
        
        fprintf('    读取文件: %s ...\n', file_path);
        try
            data = load(file_path);
            
            % 提取必要坐标轴
            if isfield(data, 'Lon'), Lon = data.Lon; elseif isfield(data, 'lon'), Lon = data.lon; else, error('无Lon'); end
            if isfield(data, 'Lat'), Lat = data.Lat; elseif isfield(data, 'lat'), Lat = data.lat; else, error('无Lat'); end
            
            if isfield(data, 'time_vec'), time_vec = data.time_vec; 
            elseif isfield(data, 'time_axis'), time_vec = data.time_axis;
            else
                warning('    未找到 time_vec，使用默认时间组 (2005-2024)');
                [Y, M] = meshgrid(2005:2024, 1:12);
                time_vec = Y(:) + (M(:)-0.5)/12; 
            end
            time_vec = time_vec(:);
            ntime = length(time_vec);
            
            % 为该态该数据开辟趋势与显著性结构体
            Trend_Detail = struct();
            Sig_Detail = struct();
            
            % 构建极速 OLS X 矩阵
            t_center = time_vec(1:ntime) - mean(time_vec(1:ntime));
            X = [t_center, ones(ntime, 1)];
            Sxx = sum(t_center.^2);
            df = ntime - 2;
            
            valid_vars_count = 0;
            
            fprintf('    开始批量极速向量化计算 %d 个分解项趋势...\n', num_vars);
            % 遍历28个变量
            for v = 1:num_vars
                v_name = var_list{v};
                if ~isfield(data, v_name), continue; end
                
                valid_vars_count = valid_vars_count + 1;
                
                % 获取 3D 数据
                D_3D = data.(v_name);
                [dim1, dim2, dim3] = size(D_3D);
                
                % 统一维度为 [Time, Lat, Lon] 便于切片
                % 假设常见的两种排列: (1) [Lon, Lat, Time]  (2) [Lat, Lon, Time]
                if dim3 == ntime
                    % 将 Time 挪到第1维: [Time, dim1, dim2]
                    D_3D = permute(D_3D, [3, 1, 2]);
                    nspace1 = dim1; nspace2 = dim2;
                elseif dim1 == ntime
                    nspace1 = dim2; nspace2 = dim3;
                else
                    warning('    [%s] 时间维度无法匹配，跳过本变量', v_name);
                    continue;
                end
                
                % 压平为 [Time, Space]
                n_space = nspace1 * nspace2;
                Y_2D = reshape(D_3D, ntime, n_space);
                
                % 剔除 NaN (只看首个点)
                valid_mask = ~isnan(Y_2D(1, :));
                
                % == 极速 OLS 计算 ==
                trend_1d = nan(1, n_space, 'single');
                sig_1d   = zeros(1, n_space, 'single');
                
                Y_valid = Y_2D(:, valid_mask);
                if ~isempty(Y_valid)
                    B = X \ Y_valid; 
                    trend_val = B(1, :);
                    
                    % 显著性
                    Y_pred = X * B;
                    Res = Y_valid - Y_pred;
                    SSE = sum(Res.^2, 1);
                    SE = sqrt((SSE / df) / Sxx);
                    
                    t_stat = abs(trend_val) ./ SE;
                    p_values = 2 * (1 - tcdf(t_stat, df));
                    
                    trend_1d(valid_mask) = single(trend_val);
                    sig_1d(valid_mask)   = single(p_values < 0.05);
                end
                
                % 还原回网格
                Trend_Detail.(v_name) = reshape(trend_1d, nspace1, nspace2);
                Sig_Detail.(v_name)   = reshape(sig_1d, nspace1, nspace2);
                
                if mod(valid_vars_count, 7) == 0
                    fprintf('      ... 已计算 %d 项 ...\n', valid_vars_count);
                end
            end
            
            fprintf('    计算完毕！有效项数: %d/28\n', valid_vars_count);
            if valid_vars_count > 0
                dataset_name = ds.name;
                save(out_file, 'Trend_Detail', 'Sig_Detail', 'Lon', 'Lat', 'time_vec', 'var_list', 'dataset_name', '-v7.3');
                fprintf('    保存至: %s\n', out_file);
            end
            
        catch ME
            fprintf('    [错误] %s\n', ME.message);
        end
    end
end
fprintf('\n========== 分解混合项趋势计算全部完成 ==========\n');
end
