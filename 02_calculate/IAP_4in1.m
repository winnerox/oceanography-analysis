%% IAP_4in1_Ultimate_Final.m
% =========================================================================
% 功能：IAP 数据集终极闭环！同步推演【平均态】与【标准参考态】
% 特性：内嵌均态计算 + 盐度经向插值 + 中点重构dz + 紧凑排版 + 规范命名
% =========================================================================
clear; clc; addpath('D:\work');

%% 1. 基础配置与网格提取 (含万能 dz)
DataRoot = 'D:\work\IAP_05_24'; TempDir = fullfile(DataRoot, 'TEMP'); SaltDir = fullfile(DataRoot, 'SALT');
OutputDir = 'D:\work\IAP_mat_data'; MeanFile = fullfile(OutputDir, 'IAP_Mean_State_05_24.mat'); 
Years = 2005:2024; Months = 1:12; MaxOrder = 8; MaxDepth = 2000;
SA_ref_val = 35.16504; CT_ref_val = 0.0; if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

fprintf('>> [1/6] 初始化 TEOS-10 引擎与提取 IAP 网格...\n'); Engine = TEOS10_General_Engine(MaxOrder);
TempSample = fullfile(TempDir, 'IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc');
SaltSample = fullfile(SaltDir, 'IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc');
Lon = mod(single(ncread(TempSample, 'lon')), 360); Lat = single(ncread(TempSample, 'lat')); Depth = single(ncread(TempSample, 'depth_std'));
SaltLon = single(ncread(SaltSample, 'lon')); NeedInterp = max(abs(Lon - mod(SaltLon,360) + 0.5)) < 0.01;
DepthIdx = find(Depth <= MaxDepth); CalcDepth = Depth(DepthIdx); num_depth = length(CalcDepth);
[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon); P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D)); [Ny, Nz, Nx] = size(P_3D);
LON_GSW = max(0.5, min(359.5, LON_3D)); LAT_GSW = max(-86, min(90, LAT_3D)); P_GSW = max(0.1, min(6131, P_3D));

% 🚨 IAP 万能 dz 重构 (中点法 + 2000m铁闸门)
try depth_bnds = double(ncread(TempSample, 'depth_bnds'));
catch
    depth_bnds = zeros(2, num_depth); depth_bnds(1, 1) = 0;
    for k = 1:num_depth-1, mid_p = (CalcDepth(k) + CalcDepth(k+1))/2; depth_bnds(2,k)=mid_p; depth_bnds(1,k+1)=mid_p; end
    depth_bnds(2, num_depth) = CalcDepth(num_depth) + (CalcDepth(num_depth) - CalcDepth(num_depth-1))/2;
end
tmp_dz = depth_bnds(2,:) - depth_bnds(1,:); cumsum_dz_1d = cumsum(tmp_dz); cumsum_dz_1d(cumsum_dz_1d > MaxDepth) = MaxDepth; 
dz_1d = cumsum_dz_1d; dz_1d(2:end) = cumsum_dz_1d(2:end) - cumsum_dz_1d(1:end-1);
dz_3D = repmat(reshape(single(dz_1d(1:num_depth)), 1, num_depth, 1), [Ny, 1, Nx]); 

%% 2. 获取或全自动计算多年平均气候态
fprintf('>> [2/6] 准备【真实平均气候态】与【标准参考态】...\n');
if exist(MeanFile, 'file'), fprintf('   ✅ 加载平均态缓存: %s\n', MeanFile); load(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean');
else
    fprintf('   ⏳ 未检测到缓存，全自动计算 IAP 20年均态...\n');
    Sum_T = zeros(Ny, Nz, Nx, 'double'); Sum_S = zeros(Ny, Nz, Nx, 'double'); Count = zeros(Ny, Nz, Nx, 'uint16');
    warning('off', 'all');
    for y = Years, for m = Months
        FileT = fullfile(TempDir, sprintf('IAPv4_Temp_monthly_1_6000m_year_%d_month_%02d.nc', y, m));
        FileS = fullfile(SaltDir, sprintf('IAPv2_Salinity_monthly_1_6000m_year_%d_month_%02d.nc', y, m));
        if ~exist(FileT, 'file') || ~exist(FileS, 'file'), continue; end
        T_raw = double(permute(ncread(FileT, 'temp', [1 1 1], [Nz Inf Inf]), [3, 1, 2]));
        S_raw = double(permute(ncread(FileS, 'salinity', [1 1 1], [Nz Inf Inf]), [3, 1, 2]));
        if NeedInterp
            S_tmp = zeros(size(T_raw), 'double'); SaltLon_ext = double([SaltLon(end)-360; SaltLon(:); SaltLon(1)+360]);
            for k=1:Nz, for j=1:Ny, S_s=squeeze(S_raw(j,k,:)); S_tmp(j,k,:)=interp1(SaltLon_ext, [S_s(end); S_s(:); S_s(1)], double(Lon), 'linear', 'extrap'); end; end
            S_raw = S_tmp;
        end
        T_raw(T_raw > 100 | T_raw < -10) = NaN; S_raw(S_raw > 100 | S_raw < 0) = NaN;
        valid = ~isnan(T_raw) & ~isnan(S_raw); Sum_T(valid) = Sum_T(valid) + T_raw(valid); Sum_S(valid) = Sum_S(valid) + S_raw(valid); Count(valid) = Count(valid) + 1;
    end; fprintf('      ... 均态已累加 %d 年\n', y); end; warning('on', 'all');
    Mean_T = single(Sum_T ./ double(Count)); Mean_S = single(Sum_S ./ double(Count)); MaskNaN_Mean = (Count == 0) | isnan(Mean_T) | isnan(Mean_S);
    save(MeanFile, 'Mean_T', 'Mean_S', 'MaskNaN_Mean', '-v7.3'); fprintf('   ✅ IAP 平均态计算完毕并缓存！\n');
end

Mean_T(MaskNaN_Mean) = 0; Mean_S(MaskNaN_Mean) = 0;
SA_Avg_3D = single(gsw_SA_from_SP(double(Mean_S), double(P_GSW), double(LON_GSW), double(LAT_GSW))); CT_Avg_3D = single(gsw_CT_from_t(SA_Avg_3D, double(Mean_T), double(P_GSW)));
SA_Std_3D = ones(Ny, Nz, Nx, 'single')*single(SA_ref_val); CT_Std_3D = ones(Ny, Nz, Nx, 'single')*single(CT_ref_val);

%% 3. 计算导数与绝对密度底座
fprintf('>> [3/6] 推演高阶泰勒系数及绝对密度底座...\n');
[Derivs_T_Avg, Derivs_S_Avg, Derivs_Cross_Avg] = Engine.calculate_all_derivatives(SA_Avg_3D, CT_Avg_3D, P_GSW, MaxOrder);
[Derivs_T_Std, Derivs_S_Std, Derivs_Cross_Std] = Engine.calculate_all_derivatives(SA_Std_3D, CT_Std_3D, P_GSW, MaxOrder);
Rho_Ref_Avg_3D = single(gsw_rho(double(SA_Avg_3D), double(CT_Avg_3D), double(P_GSW))); Rho_Ref_Std_3D = single(gsw_rho(double(SA_Std_3D), double(CT_Std_3D), double(P_GSW)));
InvFact = zeros(21,1,'single'); for i=0:20, InvFact(i+1) = single(1/factorial(i)); end

%% 4. 主循环双线推演
fprintf('>> [4/6] 启动 IAP 核心大循环...\n');
TotalSteps = length(Years)*length(Months); Time_Axis = zeros(TotalSteps, 1);
for y_idx = 1:length(Years), for m = 1:12, Time_Axis((y_idx-1)*12+m) = Years(y_idx) + (m-0.5)/12; end; end
T_Avg_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); S_Avg_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); C_Avg_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); Exact_Avg_All=nan(Ny,Nx,TotalSteps,'single');
T_Std_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); S_Std_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); C_Std_All=nan(Ny,Nx,TotalSteps,MaxOrder,'single'); Exact_Std_All=nan(Ny,Nx,TotalSteps,'single');
CrossDetail_Avg=struct(); CrossDetail_Std=struct();
for n=2:MaxOrder, for k=1:n-1, j=n-k; v=sprintf('Cross_T%dS%d',k,j); CrossDetail_Avg.(v)=nan(Ny,Nx,TotalSteps,'single'); CrossDetail_Std.(v)=nan(Ny,Nx,TotalSteps,'single'); end; end

warning('off', 'all');
for y_idx = 1:length(Years), y = Years(y_idx);
    for m = Months
        idx = (y_idx - 1) * 12 + m;
        FileT = fullfile(TempDir, sprintf('IAPv4_Temp_monthly_1_6000m_year_%d_month_%02d.nc', y, m)); FileS = fullfile(SaltDir, sprintf('IAPv2_Salinity_monthly_1_6000m_year_%d_month_%02d.nc', y, m));
        if ~exist(FileT, 'file') || ~exist(FileS, 'file'), continue; end
        try
            T_inst = single(permute(ncread(FileT, 'temp', [1 1 1], [Nz Inf Inf]), [3, 1, 2])); S_raw = single(permute(ncread(FileS, 'salinity', [1 1 1], [Nz Inf Inf]), [3, 1, 2]));
            if NeedInterp
                S_inst = zeros(size(T_inst), 'single'); SaltLon_ext = [SaltLon(end)-360; SaltLon(:); SaltLon(1)+360];
                for k=1:Nz, for j=1:Ny, S_s=squeeze(S_raw(j,k,:)); S_inst(j,k,:)=reshape(interp1(double(SaltLon_ext), double([S_s(end); S_s(:); S_s(1)]), double(Lon), 'linear', 'extrap'),1,1,[]); end; end
            else, S_inst = S_raw; end
            
            MaskInvalid = (T_inst > 100) | (T_inst < -10) | (S_inst > 100) | (S_inst < 0) | isnan(T_inst) | isnan(S_inst);
            MaskValid = ~MaskNaN_Mean & ~MaskInvalid; idx_valid = find(MaskValid); 
            HasGap_2D = squeeze(any(MaskInvalid & ~MaskNaN_Mean, 2)); Final_Mask = (squeeze(sum(MaskValid, 2)) < 1) | HasGap_2D; 
            if isempty(idx_valid), continue; end
            
            SA_1D = gsw_SA_from_SP(double(S_inst(idx_valid)), double(P_GSW(idx_valid)), double(LON_GSW(idx_valid)), double(LAT_GSW(idx_valid))); CT_1D = gsw_CT_from_t(SA_1D, double(T_inst(idx_valid)), double(P_GSW(idx_valid)));
            grav_3D = single(gsw_grav(LAT_GSW, P_GSW));
            SA_surf = gsw_SA_from_SP(double(squeeze(S_inst(:,1,:))), 0, double(squeeze(LON_GSW(:,1,:))), double(squeeze(LAT_GSW(:,1,:)))); CT_surf = gsw_CT_from_t(SA_surf, double(squeeze(T_inst(:,1,:))), 0);
            rho_a = reshape(single(gsw_rho(SA_surf, CT_surf, 0)), Ny, 1, Nx); grav_a = reshape(single(gsw_grav(double(squeeze(LAT_GSW(:,1,:))), 0)), Ny, 1, Nx);
            rho_Inst_1D = gsw_rho(SA_1D, CT_1D, double(P_GSW(idx_valid)));
            
            %% 【均态推演】
            E_3D=zeros(Ny,Nz,Nx,'single'); E_3D(idx_valid)=single(rho_Inst_1D - Rho_Ref_Avg_3D(idx_valid)); v=squeeze(-single(sum(E_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); v(Final_Mask)=NaN; Exact_Avg_All(:,:,idx)=v;
            dT_1D = CT_1D - double(CT_Avg_3D(idx_valid)); dS_1D = SA_1D - double(SA_Avg_3D(idx_valid));
            for n = 1:MaxOrder
                T_3D=zeros(Ny,Nz,Nx,'single'); T_3D(idx_valid)=single(double(InvFact(n+1)).*double(Derivs_T_Avg{n}(idx_valid)).*(dT_1D.^n));
                S_3D=zeros(Ny,Nz,Nx,'single'); S_3D(idx_valid)=single(double(InvFact(n+1)).*double(Derivs_S_Avg{n}(idx_valid)).*(dS_1D.^n)); C_1D=zeros(size(idx_valid),'double');
                if n>=2, for k=1:n-1, j=n-k; term=(double(InvFact(k+1)*InvFact(j+1))).*double(Derivs_Cross_Avg{k,j}(idx_valid)).*(dT_1D.^k).*(dS_1D.^j); C_1D=C_1D+term;
                t_3D=zeros(Ny,Nz,Nx,'single'); t_3D(idx_valid)=single(term); v=squeeze(-single(sum(t_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); v(Final_Mask)=NaN; CrossDetail_Avg.(sprintf('Cross_T%dS%d',k,j))(:,:,idx)=v; end; end
                C_3D=zeros(Ny,Nz,Nx,'single'); C_3D(idx_valid)=single(C_1D);
                vT=squeeze(-single(sum(T_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vT(Final_Mask)=NaN; T_Avg_All(:,:,idx,n)=vT;
                vS=squeeze(-single(sum(S_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vS(Final_Mask)=NaN; S_Avg_All(:,:,idx,n)=vS;
                vC=squeeze(-single(sum(C_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vC(Final_Mask)=NaN; C_Avg_All(:,:,idx,n)=vC;
            end
            
            %% 【标准态推演】
            E_3D=zeros(Ny,Nz,Nx,'single'); E_3D(idx_valid)=single(rho_Inst_1D - Rho_Ref_Std_3D(idx_valid)); v=squeeze(-single(sum(E_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); v(Final_Mask)=NaN; Exact_Std_All(:,:,idx)=v;
            dT_1D = CT_1D - double(CT_Std_3D(idx_valid)); dS_1D = SA_1D - double(SA_Std_3D(idx_valid));
            for n = 1:MaxOrder
                T_3D=zeros(Ny,Nz,Nx,'single'); T_3D(idx_valid)=single(double(InvFact(n+1)).*double(Derivs_T_Std{n}(idx_valid)).*(dT_1D.^n));
                S_3D=zeros(Ny,Nz,Nx,'single'); S_3D(idx_valid)=single(double(InvFact(n+1)).*double(Derivs_S_Std{n}(idx_valid)).*(dS_1D.^n)); C_1D=zeros(size(idx_valid),'double');
                if n>=2, for k=1:n-1, j=n-k; term=(double(InvFact(k+1)*InvFact(j+1))).*double(Derivs_Cross_Std{k,j}(idx_valid)).*(dT_1D.^k).*(dS_1D.^j); C_1D=C_1D+term;
                t_3D=zeros(Ny,Nz,Nx,'single'); t_3D(idx_valid)=single(term); v=squeeze(-single(sum(t_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); v(Final_Mask)=NaN; CrossDetail_Std.(sprintf('Cross_T%dS%d',k,j))(:,:,idx)=v; end; end
                C_3D=zeros(Ny,Nz,Nx,'single'); C_3D(idx_valid)=single(C_1D);
                vT=squeeze(-single(sum(T_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vT(Final_Mask)=NaN; T_Std_All(:,:,idx,n)=vT;
                vS=squeeze(-single(sum(S_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vS(Final_Mask)=NaN; S_Std_All(:,:,idx,n)=vS;
                vC=squeeze(-single(sum(C_3D.*grav_3D.*dz_3D,2))./(rho_a.*grav_a)*1000); vC(Final_Mask)=NaN; C_Std_All(:,:,idx,n)=vC;
            end
        catch ME, continue; end
    end; fprintf('   ... 年份 %d 数据推演完成\n', y);
end; warning('on', 'all');

%% 5. 结果保存 (国际规范命名)
fprintf('>> [5/6] 正在打包输出...\n'); lat = Lat; lon = mod(Lon,360);
SSLA_Exact_Avg = permute(Exact_Avg_All, [2, 1, 3]); TSLA_AllOrders = permute(T_Avg_All, [2, 1, 3, 4]); HSLA_AllOrders = permute(S_Avg_All, [2, 1, 3, 4]); Cross_AllOrders = permute(C_Avg_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'IAP_SSLA_Exact_Avg.mat'), 'SSLA_Exact_Avg', 'lat', 'lon', 'Time_Axis', '-v7.3'); save(fullfile(OutputDir, 'IAP_TSLA_Terms_1to8_Avg.mat'), 'TSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');
save(fullfile(OutputDir, 'IAP_HSLA_Terms_1to8_Avg.mat'), 'HSLA_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3'); save(fullfile(OutputDir, 'IAP_Cross_Terms_1to8_Avg.mat'), 'Cross_AllOrders', 'lat', 'lon', 'Time_Axis', '-v7.3');
MFinal_Avg = struct('lon', lon, 'lat', lat, 'time_vec', Time_Axis); for n=2:MaxOrder, for k=1:n-1, j=n-k; v=sprintf('Cross_T%dS%d', k, j); MFinal_Avg.(v) = permute(CrossDetail_Avg.(v), [2, 1, 3]); end; end
save(fullfile(OutputDir, 'IAP_CrossDetail_Avg.mat'), '-struct', 'MFinal_Avg', '-v7.3');

SSLA_Exact_StdRef = permute(Exact_Std_All, [2, 1, 3]); TSLA_AllOrders_Std = permute(T_Std_All, [2, 1, 3, 4]); HSLA_AllOrders_Std = permute(S_Std_All, [2, 1, 3, 4]); Cross_AllOrders_Std = permute(C_Std_All, [2, 1, 3, 4]);
save(fullfile(OutputDir, 'IAP_SSLA_Exact_StdRef.mat'), 'SSLA_Exact_StdRef', 'lat', 'lon', 'Time_Axis', '-v7.3'); save(fullfile(OutputDir, 'IAP_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3');
save(fullfile(OutputDir, 'IAP_HSLA_Terms_1to8_StdRef.mat'), 'HSLA_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3'); save(fullfile(OutputDir, 'IAP_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders_Std', 'lat', 'lon', 'Time_Axis', '-v7.3');
MFinal_Std = struct('lon', lon, 'lat', lat, 'time_vec', Time_Axis); for n=2:MaxOrder, for k=1:n-1, j=n-k; v=sprintf('Cross_T%dS%d', k, j); MFinal_Std.(v) = permute(CrossDetail_Std.(v), [2, 1, 3]); end; end
save(fullfile(OutputDir, 'IAP_CrossDetail_Std.mat'), '-struct', 'MFinal_Std', '-v7.3'); fprintf('\n>>> 🎉 IAP 终极大满贯 (含全自动均态)！ <<<\n');