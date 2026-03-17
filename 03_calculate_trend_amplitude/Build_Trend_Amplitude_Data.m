%% Build_Trend_Amplitude_Data.m
% =========================================================================
% 功能：调用 gmt_harmonic_new 批量计算三大数据集的趋势、振幅、显著性
% 目标：生成完美适配 Plot_Trend_Combined 等画图脚本的 3D 堆叠矩阵
% =========================================================================
clear; clc; close all;

%% 1. 基础配置
Datasets = {'EN4', 'IAP', 'Ishii'};
States   = {'Average', 'StdRef'}; % 这里的名字完全贴合画图脚本的要求
BaseDir  = 'D:\work';

% 创建画图脚本需要的专用结果目录
TrendDir = fullfile(BaseDir, 'Task_Convergence', 'Trend_Results');
AmpDir   = fullfile(BaseDir, 'Task_Convergence', 'Amplitude_Results');
if ~exist(TrendDir, 'dir'), mkdir(TrendDir); end
if ~exist(AmpDir, 'dir'), mkdir(AmpDir); end

fprintf('========== 🌊 启动全球调和分析与画图数据桥接器 ==========\n\n');

%% 2. 核心大循环
for d = 1:length(Datasets)
    DatasetName = Datasets{d};
    DataDir = fullfile(BaseDir, sprintf('%s_mat_data', DatasetName));
    
    for s = 1:length(States)
        StateName = States{s};
        fprintf('>> 正在处理: [%s] 数据集 -> 【%s】状态\n', DatasetName, StateName);
        
        % 匹配之前算好的文件后缀 (Average 对应 _Avg.mat, StdRef 对应 _StdRef.mat)
        if strcmp(StateName, 'Average')
            suffix = 'Avg';
        else
            suffix = 'StdRef';
        end
        
        f_SSLA = fullfile(DataDir, sprintf('%s_SSLA_Exact_%s.mat', DatasetName, suffix));
        f_TSLA = fullfile(DataDir, sprintf('%s_TSLA_Terms_1to8_%s.mat', DatasetName, suffix));
        f_HSLA = fullfile(DataDir, sprintf('%s_HSLA_Terms_1to8_%s.mat', DatasetName, suffix));
        
        if ~exist(f_SSLA, 'file') || ~exist(f_TSLA, 'file') || ~exist(f_HSLA, 'file')
            fprintf('   ⚠️ 找不到输入文件，跳过: %s %s\n', DatasetName, StateName);
            continue;
        end
        
        %% 3. 智能读取并提取数据
        fprintf('   ... 加载 20 年时空矩阵...\n');
        tmp = load(f_SSLA, 'Time_Axis', 'lat', 'lon');
        
        % 🚨 关键转换：画图脚本需要首字母大写的 Lon 和 Lat
        Lat = tmp.lat; Lon = tmp.lon; 
        
        % 🚨 关键防崩：gmt_harmonic_new 内部要求时间 t 必须是 1xN 的行向量
        Time_Vec = reshape(tmp.Time_Axis, 1, []); 
        
        SSLA_Exact_All = extract_matrix(load(f_SSLA)); % [Ny, Nx, Nt]
        TSLA_AllOrders = extract_matrix(load(f_TSLA)); % [Ny, Nx, Nt, 8]
        HSLA_AllOrders = extract_matrix(load(f_HSLA)); % [Ny, Nx, Nt, 8]
        
        [Ny, Nx, Nt, MaxOrder] = size(TSLA_AllOrders);
        
        %% 4. 初始化画图脚本所需的 3D 堆叠矩阵 [Lat, Lon, Order]
        trend_SSLA = nan(Ny, Nx, 1); sig_SSLA = nan(Ny, Nx, 1); amp_SSLA = nan(Ny, Nx, 1);
        trend_TSLA = nan(Ny, Nx, MaxOrder); sig_TSLA = nan(Ny, Nx, MaxOrder); amp_TSLA = nan(Ny, Nx, MaxOrder);
        trend_HSLA = nan(Ny, Nx, MaxOrder); sig_HSLA = nan(Ny, Nx, MaxOrder); amp_HSLA = nan(Ny, Nx, MaxOrder);
        
        fprintf('   ... 启动逐像元最小二乘调和分析 (需要运算几分钟)...\n');
        
        % --- 计算 1: 精确解 SSLA ---
        fprintf('       -> 提取 SSLA Exact 趋势与振幅...\n');
        [Amp1, ~, ~, ~, ~, ~, ~, ~, Trend, ~, Trend_sig, ~, ~, ~, ~, ~] = ...
            gmt_harmonic_new(Time_Vec, [], double(SSLA_Exact_All));
        trend_SSLA(:,:,1) = Trend; sig_SSLA(:,:,1) = Trend_sig; amp_SSLA(:,:,1) = Amp1;
        
        % --- 计算 2: T 项 (1到8阶循环提取) ---
        for n = 1:MaxOrder
            fprintf('       -> 提取 TSLA 第 %d 阶...\n', n);
            data_slice = double(TSLA_AllOrders(:,:,:,n));
            [Amp1, ~, ~, ~, ~, ~, ~, ~, Trend, ~, Trend_sig, ~, ~, ~, ~, ~] = ...
                gmt_harmonic_new(Time_Vec, [], data_slice);
            trend_TSLA(:,:,n) = Trend; sig_TSLA(:,:,n) = Trend_sig; amp_TSLA(:,:,n) = Amp1;
        end
        
        % --- 计算 3: S 项 (1到8阶循环提取) ---
        for n = 1:MaxOrder
            fprintf('       -> 提取 HSLA 第 %d 阶...\n', n);
            data_slice = double(HSLA_AllOrders(:,:,:,n));
            [Amp1, ~, ~, ~, ~, ~, ~, ~, Trend, ~, Trend_sig, ~, ~, ~, ~, ~] = ...
                gmt_harmonic_new(Time_Vec, [], data_slice);
            trend_HSLA(:,:,n) = Trend; sig_HSLA(:,:,n) = Trend_sig; amp_HSLA(:,:,n) = Amp1;
        end
        
        %% 5. 封包保存到画图脚本的指定路径
        OutTrendFile = fullfile(TrendDir, sprintf('%s_%s_Trends.mat', DatasetName, StateName));
        save(OutTrendFile, 'Lon', 'Lat', ...
            'trend_SSLA', 'sig_SSLA', ...
            'trend_TSLA', 'sig_TSLA', ...
            'trend_HSLA', 'sig_HSLA', '-v7.3');
            
        OutAmpFile = fullfile(AmpDir, sprintf('%s_%s_Amplitude.mat', DatasetName, StateName));
        save(OutAmpFile, 'Lon', 'Lat', ...
            'amp_SSLA', 'amp_TSLA', 'amp_HSLA', '-v7.3');
            
        fprintf('   ✅ [%s %s] 完美适配的画图数据已保存！\n\n', DatasetName, StateName);
    end
end

fprintf('🎉🎉 全部桥接转换完毕！现在可以直接运行 Plot_Trend_Combined 等画图脚本了！\n');

%% ======== 辅助函数：智能提取数据矩阵 ========
function DataMat = extract_matrix(StructData)
    fields = fieldnames(StructData);
    for i = 1:length(fields)
        if ndims(StructData.(fields{i})) >= 3 % 提取出 3D 或 4D 的核心矩阵
            DataMat = StructData.(fields{i});
            return;
        end
    end
    error('未在文件中找到高维数据矩阵！');
end