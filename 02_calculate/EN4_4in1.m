%% EN4_4in1_Ultimate_Diamond_Auto.m
% =========================================================================
% 功能：终极闭环！一次性并行计算 EN4 的【平均态】与【标准态】(1-8阶全家桶)
% 状态：[绝对防弹版 + 全自动均态缓存] 已排除 NaN 污染、I/O 冗余及老版本兼容性隐患
% =========================================================================
clear; clc;
addpath('D:\work');

%% 1. 基础参数与路径配置
DataRoot = 'D:\work\EN4_analyses_c13_last20years';
OutputDir = 'D:\work\EN4_mat_data'; 
MeanFile = fullfile(OutputDir, 'EN4_Mean_State.mat'); 

Years = 2005:2024;
Months = 1:12;
MaxOrder = 8;
MaxDepth = 2000;

% 标准参考态参数
SA_ref_val = 35.16504; 
CT_ref_val = 0.0;      

if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

%% 2. 初始化高阶核心引擎
fprintf('>> [1/6] 初始化 TEOS-10 升级版热力学引擎...\n');
try
    Engine = TEOS10_General_Engine(MaxOrder);
catch
    error('❌ 无法加载 TEOS10_General_Engine。');
end

%% 3. 提取原生网格与官方 depth_bnds (🚨 终极网格厚度截断)
fprintf('>> [2/6] 提取网格并生成 2000m 绝对截断积分权重...\n');
SampleFile = fullfile(DataRoot, 'EN.4.2.2.analyses.c13.2005', 'EN.4.2.2.f.analysis.c13.200501.nc');

Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, 'depth'));

DepthIdx = find(Depth <= MaxDepth); 
CalcDepth = Depth(DepthIdx);
num_depth = length(CalcDepth);

[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Ny, Nz, Nx] = size(P_3D);

LON_GSW = max(0.5, min(359.5, LON_3D));
LAT_GSW = max(-86, min(90, LAT_3D));  
P_GSW   = max(0.1, min(6131, P_3D));


depth_bnds = single(ncread(SampleFile, 'depth_bnds'));
tmp_dz = depth_bnds(2,:) - depth_bnds(1,:); 
cumsum_dz_1d = cumsum(tmp_dz);
cumsum_dz_1d(cumsum_dz_1d > MaxDepth) = MaxDepth; 
dz_1d = cumsum_dz_1d;
dz_1d(2:end) = cumsum_dz_1d(2:end) - cumsum_dz_1d(1:end-1);
% 安全扩维，防止老版本 MATLAB 报错
dz_3D = repmat(reshape(single(dz_1d(1:num_depth)), 1, num_depth, 1), [Ny, 1, Nx]); 

%% 4. Phase 1: 获取或自动计算多年平均气候态
fprintf('>> [3/6] 准备【平均态】与【标准参考态】...\n');
if exist(MeanFile, 'file')
    fprintf('   ✅ 检测到平均态缓存，直接加载: %s\n', MeanFile);
    load(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean');
else
    fprintf('   ⏳ 未检测到平均态文件，正在全自动计算 20 年平均气候态 (约需几分钟)...\n');
    Sum_T = zeros(Ny, Nz, Nx, 'double');
    Sum_S = zeros(Ny, Nz, Nx, 'double');
    Count = zeros(Ny, Nz, Nx, 'uint16');
    
    warning('off', 'all');
    for y = Years
        YearDir = fullfile(DataRoot, sprintf('EN.4.2.2.analyses.c13.%d', y));
        for m = Months
            Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
            FileList = dir(fullfile(YearDir, Pattern));
            if isempty(FileList), continue; end
            
            FileT = fullfile(FileList(1).folder, FileList(1).name);
            T_raw = double(permute(ncread(FileT, 'temperature', [1 1 1 1], [Inf Inf Nz 1]), [2, 3, 1]));
            S_raw = double(permute(ncread(FileT, 'salinity', [1 1 1 1], [Inf Inf Nz 1]), [2, 3, 1]));
            
            if nanmean(T_raw(:)) > 100, T_raw = T_raw - 273.15; end
            
            T_raw(T_raw > 100 | T_raw < -10) = NaN;
            S_raw(S_raw > 100 | S_raw < 0) = NaN;
            
            valid = ~isnan(T_raw) & ~isnan(S_raw);
            Sum_T(valid) = Sum_T(valid) + T_raw(valid);
            Sum_S(valid) = Sum_S(valid) + S_raw(valid);
            Count(valid) = Count(valid) + 1;
        end
        fprintf('      ... 平均态已累加 %d 年\n', y);
    end
    warning('on', 'all');
    
    Mean_T = single(Sum_T ./ double(Count));
    Mean_S = single(Sum_S ./ double(Count));
    MaskNaN_Mean = (Count == 0) | isnan(Mean_T) | isnan(Mean_S);
    
    save(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean', '-v7.3');
    fprintf('   ✅ 平均气候态计算完毕，已自动保存为缓存文件！\n');
end

% 构造平均态底座
Mean_T(MaskNaN_Mean) = 0; Mean_S(MaskNaN_Mean) = 0;
SA_Avg_3D = single(gsw_SA_from_SP(double(Mean_S), double(P_GSW), double(LON_GSW), double(LAT_GSW)));
CT_Avg_3D = single(gsw_CT_from_pt(SA_Avg_3D, double(Mean_T)));

% 构造标准态底座
SA_Std_3D = ones(Ny, Nz, Nx, 'single') * single(SA_ref_val);
CT_Std_3D = ones(Ny, Nz, Nx, 'single') * single(CT_ref_val);

%% 5. Phase 2: 计算两套基准的高阶导数与绝对密度场
fprintf('>> [4/6] 极速推演高阶泰勒系数及绝对密度底座...\n');
[Derivs_T_Avg, Derivs_S_Avg, Derivs_Cross_Avg] = Engine.calculate_all_derivatives(SA_Avg_3D, CT_Avg_3D, P_GSW, MaxOrder);
[Derivs_T_Std, Derivs_S_Std, Derivs_Cross_Std] = Engine.calculate_all_derivatives(SA_Std_3D, CT_Std_3D, P_GSW, MaxOrder);

% 🚨 性能护航：把基准态的 3D 绝对密度提前算出，不用每个月都去算一遍！
Rho_Ref_Avg_3D = single(gsw_rho(double(SA_Avg_3D), double(CT_Avg_3D), double(P_GSW)));
Rho_Ref_Std_3D = single(gsw_rho(double(SA_Std_3D), double(CT_Std_3D), double(P_GSW)));

InvFact = zeros(21,1,'single'); 
for i=0:20, InvFact(i+1) = single(1/factorial(i)); end

%% 6. Phase 3: 主循环 (一次读取，双线推演！)
fprintf('>> [5/6] 启动核心大循环 (1次I/O -> 2套结果同步生成)...\n');
TotalSteps = length(Years) * length(Months);
Time_Axis = zeros(TotalSteps, 1);

% ============ 🚨 安全防弹护航：全部预分配为 nan ============
T_Avg_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
S_Avg_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
C_Avg_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
Exact_Avg_All = nan(Ny, Nx, TotalSteps, 'single');
CrossDetail_Avg = struct();

T_Std_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
S_Std_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
C_Std_All = nan(Ny, Nx, TotalSteps, MaxOrder, 'single');
Exact_Std_All = nan(Ny, Nx, TotalSteps, 'single');
CrossDetail_Std = struct();

for n = 2:MaxOrder
    for k = 1:(n-1)
        j = n - k;
        var_name = sprintf('Cross_T%dS%d', k, j);
        CrossDetail_Avg.(var_name) = nan(Ny, Nx, TotalSteps, 'single');
        CrossDetail_Std.(var_name) = nan(Ny, Nx, TotalSteps, 'single');
    end
end

idx = 0; warning('off', 'all');
for y = Years
    YearDir = fullfile(DataRoot, sprintf('EN.4.2.2.analyses.c13.%d', y));
    for m = Months
        idx = idx + 1;
        Time_Axis(idx) = y + (m-0.5)/12;
        
        Pattern = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        FileList = dir(fullfile(YearDir, Pattern));
        if isempty(FileList), continue; end % 若无文件，自动保留安全的 NaN
        
        FileT = fullfile(FileList(1).folder, FileList(1).name);
        
        try
            % 提取当月温盐
            T_inst = single(permute(ncread(FileT, 'temperature', [1 1 1 1], [Inf Inf Nz 1]), [2, 3, 1]));
            S_inst = single(permute(ncread(FileT, 'salinity', [1 1 1 1], [Inf Inf Nz 1]), [2, 3, 1]));
            if mean(T_inst(:), 'omitnan') > 100, T_inst = T_inst - single(273.15); end
            
            % 严苛的 0~2000m 同源掩膜
            MaskInvalid = (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0) | isnan(T_inst) | isnan(S_inst);
            MaskValid = ~MaskNaN_Mean & ~MaskInvalid;
            idx_valid = find(MaskValid); 
            
            DynamicGap = MaskInvalid & ~MaskNaN_Mean;
            HasGap_2D = squeeze(any(DynamicGap, 2)); 
            MaskLand = (squeeze(sum(MaskValid, 2)) < 1);
            Final_Mask = (MaskLand | HasGap_2D); 
            
            if isempty(idx_valid), continue; end
            
            % 当月物理状态与边界条件 (仅算一次！)
            SA_1D = gsw_SA_from_SP(double(S_inst(idx_valid)), double(P_GSW(idx_valid)), double(LON_GSW(idx_valid)), double(LAT_GSW(idx_valid)));
            CT_1D = gsw_CT_from_pt(SA_1D, double(T_inst(idx_valid)));
            
            grav_3D = single(gsw_grav(LAT_GSW, P_GSW));
            SA_surf = gsw_SA_from_SP(double(squeeze(S_inst(:,1,:))), 0, double(squeeze(LON_GSW(:,1,:))), double(squeeze(LAT_GSW(:,1,:))));
            CT_surf = gsw_CT_from_pt(SA_surf, double(squeeze(T_inst(:,1,:))));
            rho_a  = reshape(single(gsw_rho(SA_surf, CT_surf, 0)), Ny, 1, Nx);
            grav_a = reshape(single(gsw_grav(double(squeeze(LAT_GSW(:,1,:))), 0)), Ny, 1, Nx);
            
            % 极其耗时的操作：整月实际密度
            rho_Inst_1D = gsw_rho(SA_1D, CT_1D, double(P_GSW(idx_valid)));
            
            %% ========== 【分支 A：真实平均态推演】 ==========
            rho_Ref_Avg_1D = Rho_Ref_Avg_3D(idx_valid); % 直接调用循环外算好的底座！
            E_3D = zeros(Ny, Nz, Nx, 'single'); E_3D(idx_valid) = single(rho_Inst_1D - rho_Ref_Avg_1D);
            val_E = squeeze(-single(sum(E_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
            val_E(Final_Mask) = NaN; Exact_Avg_All(:,:,idx) = val_E;
            
            dT_Avg_1D = CT_1D - double(CT_Avg_3D(idx_valid));
            dS_Avg_1D = SA_1D - double(SA_Avg_3D(idx_valid));
            
            for n = 1:MaxOrder
                T_1D = double(InvFact(n+1)) .* double(Derivs_T_Avg{n}(idx_valid)) .* (dT_Avg_1D .^ n);
                S_1D = double(InvFact(n+1)) .* double(Derivs_S_Avg{n}(idx_valid)) .* (dS_Avg_1D .^ n);
                C_1D = zeros(size(idx_valid), 'double');
                
                if n >= 2
                    for k = 1:(n-1)
                        j = n - k;
                        term_val = (double(InvFact(k+1)) * double(InvFact(j+1))) .* double(Derivs_Cross_Avg{k,j}(idx_valid)) .* (dT_Avg_1D .^ k) .* (dS_Avg_1D .^ j);
                        C_1D = C_1D + term_val;
                        
                        term_3D = zeros(Ny, Nz, Nx, 'single'); term_3D(idx_valid) = single(term_val);
                        val_term = squeeze(-single(sum(term_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                        val_term(Final_Mask) = NaN;
                        CrossDetail_Avg.(sprintf('Cross_T%dS%d', k, j))(:,:,idx) = val_term;
                    end
                end
                
                T_3D = zeros(Ny, Nz, Nx, 'single'); T_3D(idx_valid) = single(T_1D);
                S_3D = zeros(Ny, Nz, Nx, 'single'); S_3D(idx_valid) = single(S_1D);
                C_3D = zeros(Ny, Nz, Nx, 'single'); C_3D(idx_valid) = single(C_1D);
                
                val_T = squeeze(-single(sum(T_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                val_S = squeeze(-single(sum(S_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                val_C = squeeze(-single(sum(C_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                
                val_T(Final_Mask) = NaN; val_S(Final_Mask) = NaN; val_C(Final_Mask) = NaN;
                T_Avg_All(:,:,idx,n) = val_T; S_Avg_All(:,:,idx,n) = val_S; C_Avg_All(:,:,idx,n) = val_C;
            end
            
            %% ========== 【分支 B：标准参考态推演】 ==========
            rho_Ref_Std_1D = Rho_Ref_Std_3D(idx_valid); % 直接调用循环外算好的底座！
            E_3D = zeros(Ny, Nz, Nx, 'single'); E_3D(idx_valid) = single(rho_Inst_1D - rho_Ref_Std_1D);
            val_E = squeeze(-single(sum(E_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
            val_E(Final_Mask) = NaN; Exact_Std_All(:,:,idx) = val_E;
            
            dT_Std_1D = CT_1D - double(CT_Std_3D(idx_valid));
            dS_Std_1D = SA_1D - double(SA_Std_3D(idx_valid));
            
            for n = 1:MaxOrder
                T_1D = double(InvFact(n+1)) .* double(Derivs_T_Std{n}(idx_valid)) .* (dT_Std_1D .^ n);
                S_1D = double(InvFact(n+1)) .* double(Derivs_S_Std{n}(idx_valid)) .* (dS_Std_1D .^ n);
                C_1D = zeros(size(idx_valid), 'double');
                
                if n >= 2
                    for k = 1:(n-1)
                        j = n - k;
                        term_val = (double(InvFact(k+1)) * double(InvFact(j+1))) .* double(Derivs_Cross_Std{k,j}(idx_valid)) .* (dT_Std_1D .^ k) .* (dS_Std_1D .^ j);
                        C_1D = C_1D + term_val;
                        
                        term_3D = zeros(Ny, Nz, Nx, 'single'); term_3D(idx_valid) = single(term_val);
                        val_term = squeeze(-single(sum(term_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                        val_term(Final_Mask) = NaN;
                        CrossDetail_Std.(sprintf('Cross_T%dS%d', k, j))(:,:,idx) = val_term;
                    end
                end
                
                T_3D = zeros(Ny, Nz, Nx, 'single'); T_3D(idx_valid) = single(T_1D);
                S_3D = zeros(Ny, Nz, Nx, 'single'); S_3D(idx_valid) = single(S_1D);
                C_3D = zeros(Ny, Nz, Nx, 'single'); C_3D(idx_valid) = single(C_1D);
                
                val_T = squeeze(-single(sum(T_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                val_S = squeeze(-single(sum(S_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                val_C = squeeze(-single(sum(C_3D .* grav_3D .* dz_3D, 2)) ./ (rho_a.*grav_a) * 1000);
                
                val_T(Final_Mask) = NaN; val_S(Final_Mask) = NaN; val_C(Final_Mask) = NaN;
                T_Std_All(:,:,idx,n) = val_T; S_Std_All(:,:,idx,n) = val_S; C_Std_All(:,:,idx,n) = val_C;
            end
            
        catch ME
            fprintf('   ⚠️ 计算出错 %d-%02d: %s\n', y, m, ME.message);
            continue; % 出错直接跳过，保留预设的 NaN
        end
    end
    fprintf('   ... 年份 %d 数据推演完成\n', y);
end
warning('on', 'all');

%% 7. 结果保存 (规范化命名封包输出)
fprintf('>> [6/6] 所有计算完成！正在打包输出...\n');
lat = Lat; lon = Lon;

% =========================================================================
% =============== 保存【真实平均气候态】结果 ===============
% =========================================================================

% 1. 绝对精确解 (总比容 Steric -> SSLA)
SSLA_Exact_Avg  = permute(Exact_Avg_All, [2, 1, 3]); 
save(fullfile(OutputDir, 'EN4_SSLA_Exact_Avg.mat'), 'SSLA_Exact_Avg', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 2. 温度贡献项 (热容 Thermosteric -> TSLA)
TSLA_AllOrders  = permute(T_Avg_All, [2, 1, 3, 4]); 
save(fullfile(OutputDir, 'EN4_TSLA_Terms_1to8_Avg.mat'), 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 3. 盐度贡献项 (盐容 Halosteric -> HSLA)
HSLA_AllOrders  = permute(S_Avg_All, [2, 1, 3, 4]); 
save(fullfile(OutputDir, 'EN4_HSLA_Terms_1to8_Avg.mat'), 'HSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 4. 温盐混合项 (Cross)
Cross_AllOrders = permute(C_Avg_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'EN4_Cross_Terms_1to8_Avg.mat'), 'Cross_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 5. 混合项细节展开
MFinal_Avg = struct('lon', lon, 'lat', lat, 'time_vec', Time_Axis);
for n = 2:MaxOrder
    for k = 1:(n-1)
        j = n - k;
        var_name = sprintf('Cross_T%dS%d', k, j);
        MFinal_Avg.(var_name) = permute(CrossDetail_Avg.(var_name), [2, 1, 3]); 
    end
end
save(fullfile(OutputDir, 'EN4_CrossDetail_Avg.mat'), '-struct', 'MFinal_Avg', '-v7.3');
fprintf('   ✅ 真实气候态 5 个包裹已落盘 (SSLA, TSLA, HSLA 命名已严格规范)！\n');


% =========================================================================
% =============== 保存【标准参考态】结果 ===============
% =========================================================================

% 1. 绝对精确解 (总比容 Steric -> SSLA)
SSLA_Exact_StdRef = permute(Exact_Std_All, [2, 1, 3]);
save(fullfile(OutputDir, 'EN4_SSLA_Exact_StdRef.mat'), 'SSLA_Exact_StdRef', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 2. 温度贡献项 (热容 Thermosteric -> TSLA)
TSLA_AllOrders_Std    = permute(T_Std_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'EN4_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 3. 盐度贡献项 (盐容 Halosteric -> HSLA)
HSLA_AllOrders_Std    = permute(S_Std_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'EN4_HSLA_Terms_1to8_StdRef.mat'), 'HSLA_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 4. 温盐混合项 (Cross)
Cross_AllOrders_Std   = permute(C_Std_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'EN4_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3');

% 5. 混合项细节展开
MFinal_Std = struct('lon', lon, 'lat', lat, 'time_vec', Time_Axis);
for n = 2:MaxOrder
    for k = 1:(n-1)
        j = n - k;
        var_name = sprintf('Cross_T%dS%d', k, j);
        MFinal_Std.(var_name) = permute(CrossDetail_Std.(var_name), [2, 1, 3]); 
    end
end
save(fullfile(OutputDir, 'EN4_CrossDetail_Std.mat'), '-struct', 'MFinal_Std', '-v7.3');
fprintf('   ✅ 标准参考态 5 个包裹已落盘！\n');

fprintf('\n>>> 🎉 所有数据计算完毕！ <<<\n');