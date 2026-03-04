%% Calc_Exact_IAP_Std.m
% =========================================================================
% 功能：计算 IAP【标准态】下的精确热容海平面 (Exact TSLA relative to Std State)
% 物理定义：TSLA = Integral [ rho0 * (v_inst - v_std) ] dz
% 标准态：S=35, T=0
% 优化：预构建文件列表 + 减少内存操作
% =========================================================================
clear; clc;
disp('🚀 [IAP_Std] 初始化程序...');

ws = warning('off', 'all');

DataDir = 'D:\work\IAP_05_24';
OutputDir = 'D:\work\IAP_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; rho0 = 1035.0; MaxDepth = 2000;
Std_S_Val = 35.0; Std_T_Val = 0.0;

List = dir(fullfile(DataDir, 'TEMP', '*.nc'));
if isempty(List), error('找不到 IAP 数据文件'); end
SampleFile = fullfile(List(1).folder, List(1).name);

Lon = double(ncread(SampleFile, 'lon')); 
Lat = double(ncread(SampleFile, 'lat'));
vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth_std')), vDepth='depth_std';
elseif any(strcmp(vnames, 'depth')), vDepth='depth'; 
else, vDepth='depth_std'; end
Depth = double(ncread(SampleFile, vDepth)); 

Lon(Lon > 180) = Lon(Lon > 180) - 360;
Lat(Lat < -85) = -85; Lat(Lat > 85) = 85;

DepthIdx = find(Depth <= MaxDepth); 
CalcDepth = Depth(DepthIdx);
Nx = length(Lon); Ny = length(Lat); Nz = length(DepthIdx);

[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);

dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_3D = reshape(dz, 1, 1, Nz);

disp('计算标准态密度...');
S_Std_1D = repmat(Std_S_Val, Nx*Ny*Nz, 1);
T_Std_1D = repmat(Std_T_Val, Nx*Ny*Nz, 1);
SA_Std_1D = gsw_SA_from_SP(S_Std_1D, P_3D(:), LON_3D(:), LAT_3D(:));
CT_Std_1D = gsw_CT_from_pt(SA_Std_1D, T_Std_1D);
v_Std_1D = gsw_specvol(SA_Std_1D, CT_Std_1D, P_3D(:));
rho_Std_1D = 1 ./ v_Std_1D;

TotalSteps = length(Years) * 12;
TSLA_Exact_Std = zeros(Nx, Ny, TotalSteps, 'single');
time_vec = zeros(TotalSteps, 1);

disp('预构建文件列表...');
FileList_T = cell(TotalSteps, 1);
FileList_S = cell(TotalSteps, 1);
idx = 0;
for y = Years
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        Pat = sprintf('*year_%d_month_%02d.nc', y, m);
        sT = dir(fullfile(DataDir, 'TEMP', Pat));
        sS = dir(fullfile(DataDir, 'SALT', Pat));
        if ~isempty(sT), FileList_T{idx} = fullfile(sT(1).folder, sT(1).name); end
        if ~isempty(sS), FileList_S{idx} = fullfile(sS(1).folder, sS(1).name); end
    end
end

disp('开始处理逐月时间序列...');
for idx = 1:TotalSteps
    y = floor(time_vec(idx));
    m = round((time_vec(idx) - y) * 12 + 0.5);
    fprintf('  -> %04d-%02d: ', y, m);
    
    if isempty(FileList_T{idx}) || isempty(FileList_S{idx})
        TSLA_Exact_Std(:,:,idx) = NaN;
        fprintf('文件缺失\n');
        continue;
    end
    
    try
        T_raw = double(permute(ncread(FileList_T{idx}, 'temp', [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
        S_raw = double(permute(ncread(FileList_S{idx}, 'salinity', [1 1 1], [Nz Inf Inf]), [2, 3, 1]));
        
        mask_inst = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
        idx_inst = find(mask_inst);
        
        if isempty(idx_inst)
            TSLA_Exact_Std(:,:,idx) = NaN;
            fprintf('无有效数据\n');
            continue;
        end
        
        SA_Inst = gsw_SA_from_SP(S_raw(idx_inst), P_3D(idx_inst), LON_3D(idx_inst), LAT_3D(idx_inst));
        CT_Inst = gsw_CT_from_pt(SA_Inst, T_raw(idx_inst));
        v_Inst = gsw_specvol(SA_Inst, CT_Inst, P_3D(idx_inst));
        rho_Inst = 1 ./ v_Inst;
        
        d_rho = zeros(Nx, Ny, Nz);
        d_rho(idx_inst) = rho_Inst - rho_Std_1D(idx_inst);
        d_rho(~mask_inst) = NaN;
        
        val_mm = -nansum(d_rho .* dz_3D, 3) / rho0 * 1000;
        val_mm(all(isnan(d_rho), 3)) = NaN;
        
        TSLA_Exact_Std(:,:,idx) = single(val_mm);
        fprintf('完成\n');
    catch ME
        fprintf('错误: %s\n', ME.message);
        TSLA_Exact_Std(:,:,idx) = NaN;
    end
end

SaveName = fullfile(OutputDir, 'IAP_Formula11_Exact_Std.mat');
lon = Lon; lat = Lat;
save(SaveName, 'TSLA_Exact_Std', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('已保存: %s\n', SaveName);
warning(ws);
