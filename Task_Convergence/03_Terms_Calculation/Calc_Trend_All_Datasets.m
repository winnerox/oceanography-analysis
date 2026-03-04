%% Calc_Trend_All_Datasets.m
% =========================================================================
% 计算三个数据集（EN4, IAP, Ishii）的趋势（Trend）
% 采用普通最小二乘法 (OLS) 对每个格点的时间序列进行线性回归
% 处理状态：平均态 (Average)、标准态 (StdRef)
% 处理分量：TSLA, SSLA, Cross 的 1到8阶
% =========================================================================
clear; clc;

BaseDirs = {'D:\work\EN4_TSLA_Terms', 'D:\work\IAP_mat_data', 'D:\work\Ishii_mat_data'};
Models = {'EN4', 'IAP', 'Ishii'};
States = {'Average', 'StdRef'};
MaxOrder = 8;

for m = 1:length(Models)
    Model = Models{m};
    Folder = BaseDirs{m};
    
    for s = 1:length(States)
        State = States{s};
        fprintf('\n>>> 正在处理 %s - %s 态...\n', Model, State);
        
        file_tsla  = fullfile(Folder, sprintf('%s_TSLA_Terms_1to8_%s.mat', Model, State));
        file_ssla  = fullfile(Folder, sprintf('%s_S_Terms_1to8_%s.mat', Model, State));
        file_cross = fullfile(Folder, sprintf('%s_Cross_Terms_1to8_%s.mat', Model, State));
        
        if ~exist(file_tsla, 'file')
            fprintf('  [警告] 找不到文件 %s, 跳过当前处理。\n', file_tsla);
            continue; 
        end
        
        % 加载网格坐标和时间轴
        fprintf('  加载坐标和时间轴...\n');
        info = load(file_tsla, 'Time_Axis', 'time_axis', 'time_vec', 'Lon', 'lon', 'Lat', 'lat');
        if isfield(info, 'Time_Axis'), t_vec = info.Time_Axis(:);
        elseif isfield(info, 'time_axis'), t_vec = info.time_axis(:);
        elseif isfield(info, 'time_vec'), t_vec = info.time_vec(:);
        else, error('未找到时间轴变量！');
        end
        
        if isfield(info, 'Lon'), Lon = info.Lon; else, Lon = info.lon; end
        if isfield(info, 'Lat'), Lat = info.Lat; else, Lat = info.lat; end
        
        % 1. 计算 TSLA Trend
        fprintf('  计算 TSLA 趋势 (OLS)...\n');
        load(file_tsla, 'TSLA_AllOrders');
        TSLA_Trend = calc_ols_trend(TSLA_AllOrders, t_vec, MaxOrder);
        clear TSLA_AllOrders
        
        % 2. 计算 SSLA Trend
        fprintf('  计算 SSLA 趋势 (OLS)...\n');
        if exist(file_ssla, 'file')
            load(file_ssla, 'SSLA_AllOrders');
            SSLA_Trend = calc_ols_trend(SSLA_AllOrders, t_vec, MaxOrder);
            clear SSLA_AllOrders
        else
            SSLA_Trend = [];
        end
        
        % 3. 计算 Cross Trend
        fprintf('  计算 Cross 趋势 (OLS)...\n');
        if exist(file_cross, 'file')
            load(file_cross, 'Cross_AllOrders');
            Cross_Trend = calc_ols_trend(Cross_AllOrders, t_vec, MaxOrder);
            clear Cross_AllOrders
        else
            Cross_Trend = [];
        end
        
        % 保存结果
        out_file = fullfile(Folder, sprintf('%s_Trend_1to8_%s.mat', Model, State));
        fprintf('  保存趋势结果至: %s\n', out_file);
        save(out_file, 'TSLA_Trend', 'SSLA_Trend', 'Cross_Trend', 'Lon', 'Lat', '-v7.3');
    end
end
fprintf('\n>>> 🎉 所有数据集均已完成最小二乘法(OLS)趋势计算并保存！\n');

%% Helper Function: OLS Vectorized-Loop
function TrendMap = calc_ols_trend(data4d, t_vec, max_order)
    % data4d 的预期维度为 [Lon, Lat, Time, Order]
    dim = size(data4d);
    if length(dim) == 3
        dim(4) = 1; % 可能只有 1 阶的情形
    end
    Nx = dim(1); Ny = dim(2); Nt = dim(3); MaxO = dim(4);
    
    if MaxO < max_order, max_order = MaxO; end
    
    TrendMap = nan(Nx, Ny, max_order);
    
    for k = 1:max_order
        % 将空间维度展平 [Nx*Ny, Nt]
        Y = reshape(data4d(:,:,:,k), Nx*Ny, Nt);
        
        % 找到包含有效数据（非NaN且长度达标）的格点
        valid_mask = sum(~isnan(Y), 2) > Nt * 0.6;
        valid_idx = find(valid_mask);
        
        trend_flat = nan(Nx*Ny, 1);
        
        % 中心化时间轴以解决矩阵条件数过大的警告 (ill-conditioned matrix)
        t_mean = nanmean(t_vec);
        t_centered = t_vec - t_mean;
        
        % 循环遍历所有的有效网格点以计算趋势
        % 由于只对有效点进行 OLS() 操作，MATLAB 的 for 循环在这个规模（几万点）耗时可控（几秒内）
        for i = 1:length(valid_idx)
            idx = valid_idx(i);
            y = Y(idx, :)';
            not_nan = ~isnan(y);
            
            % 仅仅对非 NaN 数据点做线性回归
            if sum(not_nan) > Nt * 0.6
                % 构造自变量矩阵 (截距列 + 中心化时间列)
                x_sub = [ones(sum(not_nan), 1), t_centered(not_nan)];
                
                % OLS 回归公式求解系数 beta = (X'*X) \ X'*y
                b = (x_sub'*x_sub) \ (x_sub'*y(not_nan));
                
                % 提取斜率 (也就是Trend)
                trend_flat(idx) = b(2);
            end
        end
        % 将一维结果重塑回 [Lon, Lat] 结构
        TrendMap(:,:,k) = reshape(trend_flat, Nx, Ny);
    end
end
