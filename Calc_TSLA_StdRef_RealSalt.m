%% Calc_TSLA_StdRef_RealSalt.m
% =========================================================================
% 标准态重算版 (串行安全版)
% 
% 核心逻辑：
%   1. S_ref = S_real (真实盐度)
%   2. T_ref = 0 (标准温度)
%   3. 每一时刻重新计算导数 (耗时操作)
%   4. 【改进】支持断点续传：每算完一年自动保存，防止崩溃丢失数据。
% =========================================================================
clear; clc;

%% 1. 配置
DataDir = 'D:\work\EN4_analyses_c13_last20years'; 
OutputDir = 'D:\work\EN4_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; 
Months = 1:12;
MaxOrder = 8;        
MaxDepth = 2000;     
rho0 = 1035.0;       
Std_Temp_Val = 0.0;  % 参考态 T=0

FinalSaveName = fullfile(OutputDir, 'EN4_TSLA_Terms_1to8_StdRef_RealSalt.mat');

%% 2. 初始化
fprintf('>> [1/4] 初始化引擎与网格...\n');
try
    Engine = TEOS10_General_Engine(); 
catch
    error('❌ 未找到 TEOS10_General_Engine，请检查路径。');
end

% 自动搜索网格文件
SampleFile = '';
for y = Years
    FList = dir(fullfile(DataDir, sprintf('*%d*',y), '*.nc'));
    if ~isempty(FList)
        SampleFile = fullfile(FList(1).folder, FList(1).name); break; 
    end
end
if isempty(SampleFile), error('❌ 未找到数据文件'); end

Lon = ncread(SampleFile, 'lon');
Lat = ncread(SampleFile, 'lat');
Depth = ncread(SampleFile, 'depth'); 
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

% 垂直权重
dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

% 预计算阶乘倒数
factorial_inv = zeros(1, MaxOrder);
for n = 1:MaxOrder, factorial_inv(n) = 1/factorial(n); end

%% 3. 遍历计算 (分年保存模式)
fprintf('>> [2/4] 开始串行遍历 (Result will be saved year by year)...\n');

for y = Years
    fprintf('>> 正在处理年份: %d ... ', y);
    TicYear = tic;
    
    % 定义该年的临时文件名
    YearFile = fullfile(OutputDir, sprintf('Temp_TSLA_%d.mat', y));
    
    % 如果这个文件已经存在，说明之前跑过了，直接跳过
    if exist(YearFile, 'file')
        fprintf('[已存在，跳过]\n');
        continue;
    end
    
    % 初始化该年的数据容器 [Nx, Ny, 12, MaxOrder]
    Year_TSLA = nan(Nx, Ny, 12, MaxOrder, 'single'); % 使用 single 节省一半内存
    
    FPath = dir(fullfile(DataDir, sprintf('*%d*',y), '*.nc'));
    if isempty(FPath), fprintf('[无数据]\n'); continue; end
    
    % 月份循环
    for m = 1:12
        % 找对应月份的文件
        FilePattern = sprintf('*%d%02d*.nc', y, m);
        TargetList = dir(fullfile(FPath(1).folder, FilePattern));
        if isempty(TargetList), continue; end
        FullPath = fullfile(TargetList(1).folder, TargetList(1).name);
        
        try
            T_3D = ncread(FullPath, 'temperature');
            S_3D = ncread(FullPath, 'salinity');
            T_3D = double(T_3D(:,:,DepthIdx)); % 转回 double 计算精度
            S_3D = double(S_3D(:,:,DepthIdx));
            if nanmean(T_3D(:)) > 100, T_3D = T_3D - 273.15; end
        catch
            continue;
        end
        
        ValidMask = ~isnan(T_3D) & ~isnan(S_3D);
        valid_idx = find(ValidMask);
        
        if ~isempty(valid_idx)
            % 提取向量
            P_vec = P_3D(valid_idx);
            S_vec = S_3D(valid_idx);
            T_vec = T_3D(valid_idx);
            
            % 1. 坐标转换
            SA_vec = gsw_SA_from_SP(S_vec, P_vec, LON_3D(valid_idx), LAT_3D(valid_idx));
            
            % 2. 构造参考态 (S=Real, T=0)
            Pt_Ref_vec = Std_Temp_Val * ones(size(SA_vec));
            CT_Ref_vec = gsw_CT_from_pt(SA_vec, Pt_Ref_vec);
            
            % 3. 构造真实态 (S=Real, T=Real)
            CT_Inst_vec = gsw_CT_from_pt(SA_vec, T_vec);
            
            % 4. 计算导数 (最耗时步骤)
            % 因为是串行，这里直接调用
            Derivs = Engine.calculate_all_orders(SA_vec, CT_Ref_vec, P_vec, MaxOrder);
            
            % 5. 泰勒求和
            Delta_Theta = CT_Inst_vec - CT_Ref_vec;
            Delta_Pow   = Delta_Theta; % 一次方
            
            for n = 1:MaxOrder
                FieldName = sprintf('d%d_T', n);
                d_rho_n = factorial_inv(n) .* Derivs.(FieldName) .* Delta_Pow;
                
                % 填回网格并积分
                Grid_Temp = zeros(Nx, Ny, Nz);
                Grid_Temp(valid_idx) = d_rho_n;
                
                % 垂直积分 TSLA (mm)
                val_mm = -(sum(Grid_Temp .* dz_perm, 3) / rho0) * 1000;
                
                Year_TSLA(:,:,m,n) = single(val_mm); % 存为 single
                
                if n < MaxOrder
                    Delta_Pow = Delta_Pow .* Delta_Theta; % 准备下一阶
                end
            end
        end
    end
    
    % 保存该年的数据 (Checkpoint)
    save(YearFile, 'Year_TSLA', 'y');
    fprintf('完成 (耗时 %.1f 分钟)\n', toc(TicYear)/60);
end

%% 4. 合并数据
fprintf('>> [3/4] 所有年份计算完毕，开始合并数据...\n');

TotalSteps = length(Years) * 12;
TSLA_AllOrders = zeros(Nx, Ny, TotalSteps, MaxOrder, 'single');
time_axis = zeros(TotalSteps, 1);

idx = 0;
for y = Years
    YearFile = fullfile(OutputDir, sprintf('Temp_TSLA_%d.mat', y));
    if exist(YearFile, 'file')
        D = load(YearFile); % 加载 Year_TSLA
        for m = 1:12
            idx = idx + 1;
            time_axis(idx) = y + (m-0.5)/12;
            TSLA_AllOrders(:,:,idx,:) = D.Year_TSLA(:,:,m,:);
        end
    else
        % 如果某一年缺失，填充 NaN
        for m = 1:12
            idx = idx + 1;
            time_axis(idx) = y + (m-0.5)/12;
            TSLA_AllOrders(:,:,idx,:) = NaN;
        end
    end
end

%% 5. 最终保存
fprintf('>> [4/4] 保存最终合并文件...\n');
lon = Lon; lat = Lat;
save(FinalSaveName, 'TSLA_AllOrders', 'lon', 'lat', 'time_axis', '-v7.3');

% 清理临时文件 (可选，建议确认无误后再手动删)
% delete(fullfile(OutputDir, 'Temp_TSLA_*.mat'));

fprintf('>> ✅ 任务全部完成！\n');
fprintf('   结果保存在: %s\n', FinalSaveName);