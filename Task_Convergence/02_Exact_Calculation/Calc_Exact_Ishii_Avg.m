%% Calc_Exact_Ishii_Avg.m
% =========================================================================
% 功能：计算 Ishii【平均态】下的精确热容海平面 (Exact TSLA relative to Mean State)
% 物理定义：TSLA = Integral [ rho0 * (v_inst - v_mean) ] dz
% 优化：预构建文件列表 + 减少内存操作
% =========================================================================
clear; clc;
disp('[Ishii_Avg] 初始化程序...');

ws = warning('off', 'all');

DataDir = 'D:\work\Ishii_05_24';
OutputDir = 'D:\work\Ishii_TSLA_Terms';
MeanFile = fullfile(OutputDir, 'Ishii_Mean_State.mat');
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; rho0 = 1035.0; MaxDepth = 2000;

if ~exist(MeanFile, 'file'), error('请先运行 Calc_Mean_Ishii.m'); end
load(MeanFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');

Lon = double(Lon); Lat = double(Lat);
if any(isnan(Lon)), Lon(isnan(Lon)) = 0; end
if any(isnan(Lat)), Lat(isnan(Lat)) = 0; end
Lat(Lat < -85) = -85; Lat(Lat > 85) = 85;

CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_3D = reshape(dz, 1, 1, Nz);

disp('计算平均态密度...');
mask_mean = (Mean_T > -100) & (Mean_S >= 0) & (Mean_S < 100);
idx_mean = find(mask_mean);

SA_Mean_1D = gsw_SA_from_SP(Mean_S(idx_mean), P_3D(idx_mean), LON_3D(idx_mean), LAT_3D(idx_mean));
CT_Mean_1D = gsw_CT_from_pt(SA_Mean_1D, Mean_T(idx_mean));
v_Mean_1D = gsw_specvol(SA_Mean_1D, CT_Mean_1D, P_3D(idx_mean));
rho_Mean_1D = 1 ./ v_Mean_1D;

rho_Mean = NaN(Nx, Ny, Nz);
rho_Mean(idx_mean) = rho_Mean_1D;

TotalSteps = length(Years) * 12;
TSLA_Exact_Avg = zeros(Nx, Ny, TotalSteps, 'single');
time_vec = zeros(TotalSteps, 1);

disp('预构建文件列表...');
FileList_T = cell(TotalSteps, 1);
FileList_S = cell(TotalSteps, 1);
month_idx = zeros(TotalSteps, 1);
idx = 0;
for y = Years
    FT = fullfile(DataDir, 'Temperature', sprintf('temp.%d.nc', y));
    FS = fullfile(DataDir, 'Salinity', sprintf('sal.%d.nc', y));
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        month_idx(idx) = m;
        if exist(FT, 'file'), FileList_T{idx} = FT; end
        if exist(FS, 'file'), FileList_S{idx} = FS; end
    end
end

disp('开始处理逐月时间序列...');
for idx = 1:TotalSteps
    y = floor(time_vec(idx));
    m = month_idx(idx);
    fprintf('  -> %04d-%02d: ', y, m);
    
    if isempty(FileList_T{idx}) || isempty(FileList_S{idx})
        TSLA_Exact_Avg(:,:,idx) = NaN;
        fprintf('文件缺失\n');
        continue;
    end
    
    try
        T_raw = double(ncread(FileList_T{idx}, 'temp', [1 1 1 m], [Inf Inf Nz 1]));
        S_raw = double(ncread(FileList_S{idx}, 'sal', [1 1 1 m], [Inf Inf Nz 1]));
        
        mask_inst = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
        idx_inst = find(mask_inst);
        
        if isempty(idx_inst)
            TSLA_Exact_Avg(:,:,idx) = NaN;
            fprintf('无有效数据\n');
            continue;
        end
        
        SA_Inst = gsw_SA_from_SP(S_raw(idx_inst), P_3D(idx_inst), LON_3D(idx_inst), LAT_3D(idx_inst));
        CT_Inst = gsw_CT_from_pt(SA_Inst, T_raw(idx_inst));
        v_Inst = gsw_specvol(SA_Inst, CT_Inst, P_3D(idx_inst));
        rho_Inst = 1 ./ v_Inst;
        
        d_rho = zeros(Nx, Ny, Nz);
        d_rho(idx_inst) = rho_Inst - rho_Mean(idx_inst);
        d_rho(~mask_mean) = NaN;
        
        val_mm = -nansum(d_rho .* dz_3D, 3) / rho0 * 1000;
        val_mm(all(isnan(d_rho), 3)) = NaN;
        
        TSLA_Exact_Avg(:,:,idx) = single(val_mm);
        fprintf('完成\n');
    catch ME
        fprintf('错误: %s\n', ME.message);
        TSLA_Exact_Avg(:,:,idx) = NaN;
    end
end

SaveName = fullfile(OutputDir, 'Ishii_Formula11_Exact_Avg.mat');
lon = Lon; lat = Lat;
save(SaveName, 'TSLA_Exact_Avg', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('已保存: %s\n', SaveName);
warning(ws);
