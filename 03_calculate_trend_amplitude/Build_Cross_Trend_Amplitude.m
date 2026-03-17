%% Build_Cross_Trend_Amplitude.m
% =========================================================================
% 功能：调用 gmt_harmonic_new 批量计算三大数据集28个交叉混合项的趋势、振幅、显著性
% 目标：生成适配 Plot_CrossDetail 等画图脚本的 3D 堆叠矩阵 [Lat, Lon, 28]
% =========================================================================
clear; clc; close all;

%% 1. 基础配置
Datasets = {'EN4', 'IAP', 'Ishii'};
States   = {'Average', 'StdRef'}; 
BaseDir  = 'D:\work';

% 创建结果目录
TrendDir = fullfile(BaseDir, 'Task_Convergence', 'Trend_Results');
AmpDir   = fullfile(BaseDir, 'Task_Convergence', 'Amplitude_Results');
if ~exist(TrendDir, 'dir'), mkdir(TrendDir); end
if ~exist(AmpDir, 'dir'), mkdir(AmpDir); end

% 28个混合项的严谨排序
cross_terms = {
    'Cross_T1S1'; 'Cross_T1S2'; 'Cross_T1S3'; 'Cross_T1S4'; 'Cross_T1S5'; 'Cross_T1S6'; 'Cross_T1S7';
    'Cross_T2S1'; 'Cross_T2S2'; 'Cross_T2S3'; 'Cross_T2S4'; 'Cross_T2S5'; 'Cross_T2S6';
    'Cross_T3S1'; 'Cross_T3S2'; 'Cross_T3S3'; 'Cross_T3S4'; 'Cross_T3S5';
    'Cross_T4S1'; 'Cross_T4S2'; 'Cross_T4S3'; 'Cross_T4S4';
    'Cross_T5S1'; 'Cross_T5S2'; 'Cross_T5S3';
    'Cross_T6S1'; 'Cross_T6S2';
    'Cross_T7S1'
};

fprintf('========== 🌊 启动混合项 (Cross) 调和分析桥接器 ==========\n');
fprintf('⚠️ 提示：需提取 28x3x2 = 168 个矩阵，这可能需要喝杯咖啡的时间...\n\n');

%% 2. 核心大循环
for d = 1:length(Datasets)
    DatasetName = Datasets{d};
    DataDir = fullfile(BaseDir, sprintf('%s_mat_data', DatasetName));
    
    for s = 1:length(States)
        StateName = States{s};
        fprintf('>> 正在处理: [%s] 数据集 -> 【%s】状态\n', DatasetName, StateName);
        
        if strcmp(StateName, 'Average'), suffix = 'Avg'; else, suffix = 'Std'; end
        f_in = fullfile(DataDir, sprintf('%s_CrossDetail_%s.mat', DatasetName, suffix));
        
        if ~exist(f_in, 'file')
            fprintf('   ⚠️ 找不到输入文件，跳过: %s\n', f_in); continue;
        end
        
        %% 3. 加载并初始化
        fprintf('   ... 加载基础数据...\n');
        tmp = load(f_in);
        
        % 规范化经纬度和时间
        Lon = tmp.lon; Lat = tmp.lat; 
        Time_Vec = reshape(tmp.time_vec, 1, []); % 必须是 1xN 行向量
        
        nx = length(Lon); ny = length(Lat);
        trend_Cross = nan(ny, nx, 28);
        sig_Cross   = nan(ny, nx, 28);
        amp_Cross   = nan(ny, nx, 28);
        
        %% 4. 逐个提取 28 个交叉项
        for i = 1:28
            term_name = cross_terms{i};
            if mod(i, 7) == 1
                fprintf('       -> 正在提取 %s (进度: %d/28)...\n', term_name, i);
            end
            
            % 读取并转置为 [Lat, Lon, Time] 适应 gmt_harmonic_new
            raw_data = tmp.(term_name); 
            data_slice = permute(raw_data, [2, 1, 3]); 
            
            % 调用老师的调和分析引擎
            [Amp1, ~, ~, ~, ~, ~, ~, ~, Trend, ~, Trend_sig, ~, ~, ~, ~, ~] = ...
                gmt_harmonic_new(Time_Vec, [], double(data_slice));
                
            trend_Cross(:,:,i) = Trend;
            sig_Cross(:,:,i)   = Trend_sig;
            amp_Cross(:,:,i)   = Amp1;
        end
        
        %% 5. 封包保存
        OutTrendFile = fullfile(TrendDir, sprintf('%s_%s_Cross_Trends.mat', DatasetName, StateName));
        save(OutTrendFile, 'Lon', 'Lat', 'trend_Cross', 'sig_Cross', 'cross_terms', '-v7.3');
            
        OutAmpFile = fullfile(AmpDir, sprintf('%s_%s_Cross_Amplitude.mat', DatasetName, StateName));
        save(OutAmpFile, 'Lon', 'Lat', 'amp_Cross', 'cross_terms', '-v7.3');
            
        fprintf('   ✅ [%s %s] Cross 画图专属矩阵已打包保存！\n\n', DatasetName, StateName);
    end
end

fprintf('🎉🎉 全部桥接转换完毕！现在可以实现秒级画图了！\n');