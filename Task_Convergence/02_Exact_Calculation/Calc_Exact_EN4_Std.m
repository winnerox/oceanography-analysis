%% Calc_Exact_EN4_Std.m
% =========================================================================
% 功能：计算 EN4【标准态】下的精确热容海平面 (Exact TSLA relative to Std State)
% 物理定义：TSLA = Integral [ rho0 * (v_inst - v_std) ] dz
% 标准态：S=35, T=0
% =========================================================================
clear; clc;
disp('🚀 [EN4_Std] 初始化程序...');

ws = warning('off', 'all');

DataDir = 'D:\work\EN4_analyses_c13_last20years';
OutputDir = 'D:\work\EN4_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; rho0 = 1035.0; MaxDepth = 2000;
Std_S_Val = 35.0; Std_T_Val = 0.0;

YearDir = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', Years(1)));
List = dir(fullfile(YearDir, '*.nc'));
if isempty(List), error('❌ 找不到 EN4 数据文件: %s', YearDir); end
SampleFile = fullfile(List(1).folder, List(1).name);

Lon = double(ncread(SampleFile, 'lon')); 
Lat = double(ncread(SampleFile, 'lat'));
vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth')), vDepth='depth'; 
elseif any(strcmp(vnames, 'deptht')), vDepth='deptht'; 
else, vDepth='depth'; end
Depth = double(ncread(SampleFile, vDepth)); 

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
dz_perm = repmat(reshape(dz, 1, 1, Nz), [Nx, Ny, 1]);

disp('⏳ 计算标准态密度...');
S_Std_1D = repmat(Std_S_Val, Nx*Ny*Nz, 1);
T_Std_1D = repmat(Std_T_Val, Nx*Ny*Nz, 1);
SA_Std_1D = gsw_SA_from_SP(S_Std_1D, P_3D(:), LON_3D(:), LAT_3D(:));
CT_Std_1D = gsw_CT_from_pt(SA_Std_1D, T_Std_1D);
v_Std_1D = gsw_specvol(SA_Std_1D, CT_Std_1D, P_3D(:));
rho_Std_1D = 1 ./ v_Std_1D;
rho_Std = reshape(rho_Std_1D, Nx, Ny, Nz);

TotalSteps = length(Years) * 12;
TSLA_Exact_Std = zeros(Nx, Ny, TotalSteps, 'single');
time_vec = zeros(TotalSteps, 1);
idx = 0;

disp('🌊 开始处理逐月时间序列...');
for y = Years
    YearDir = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y));
    
    if ~exist(YearDir, 'dir')
        fprintf('  [跳过] 找不到 %d 年目录: %s\n', y, YearDir);
        for m = 1:12
            idx = idx + 1;
            time_vec(idx) = y + (m-0.5)/12;
            TSLA_Exact_Std(:,:,idx) = NaN;
        end
        continue;
    end
    
    for m = 1:12
        idx = idx + 1;
        time_vec(idx) = y + (m-0.5)/12;
        fprintf('  -> %04d-%02d: ', y, m);
        
        FilePat = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        TargetFile = dir(fullfile(YearDir, FilePat));
        if isempty(TargetFile)
            TSLA_Exact_Std(:,:,idx) = NaN;
            fprintf('文件缺失\n');
            continue;
        end
        FullPath = fullfile(TargetFile(1).folder, TargetFile(1).name);
        
        try
            T_raw = double(ncread(FullPath, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
            S_raw = double(ncread(FullPath, 'salinity', [1 1 1 1], [Inf Inf Nz 1]));
            
            if nanmean(T_raw(:)) > 100, T_raw = T_raw - 273.15; end
            
            mask_inst = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
            idx_inst = find(mask_inst);
            
            SA_Inst_1D = gsw_SA_from_SP(S_raw(idx_inst), P_3D(idx_inst), LON_3D(idx_inst), LAT_3D(idx_inst));
            CT_Inst_1D = gsw_CT_from_pt(SA_Inst_1D, T_raw(idx_inst));
            v_Inst_1D = gsw_specvol(SA_Inst_1D, CT_Inst_1D, P_3D(idx_inst));
            rho_Inst_1D = 1 ./ v_Inst_1D;
            
            rho_Inst = NaN(Nx, Ny, Nz);
            rho_Inst(idx_inst) = rho_Inst_1D;
            
            d_rho = rho_Inst - rho_Std;
            val_mm = -nansum(d_rho .* dz_perm, 3) / rho0 * 1000;
            val_mm(all(isnan(d_rho), 3)) = NaN;
            
            TSLA_Exact_Std(:,:,idx) = single(val_mm);
            fprintf('完成\n');
        catch ME
            fprintf('错误: %s\n', ME.message);
            TSLA_Exact_Std(:,:,idx) = NaN;
        end
    end
end

SaveName = fullfile(OutputDir, 'EN4_Formula11_Exact_Std.mat');
lon = Lon; lat = Lat;
save(SaveName, 'TSLA_Exact_Std', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('🎉 已保存: %s\n', SaveName);
warning(ws);
