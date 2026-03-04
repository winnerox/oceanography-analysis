%% Calc_Mean_Ishii.m
% =========================================================================
% 功能：计算 Ishii 的平均态
% 输出：D:\work\Ishii_TSLA_Terms\Ishii_Mean_State.mat
% =========================================================================
clear; clc;

ws = warning('off', 'all');

DataDir = 'D:\work\Ishii_05_24';
OutputDir = 'D:\work\Ishii_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; MaxDepth = 2000;

List = dir(fullfile(DataDir, 'Temperature', '*.nc'));
SampleFile = fullfile(List(1).folder, List(1).name);

vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'lon')), vLon='lon'; elseif any(strcmp(vnames, 'long')), vLon='long'; elseif any(strcmp(vnames, 'longitude')), vLon='longitude'; else, vLon=vnames{1}; end
if any(strcmp(vnames, 'lat')), vLat='lat'; elseif any(strcmp(vnames, 'latitude')), vLat='latitude'; else, vLat=vnames{2}; end
if any(strcmp(vnames, 'depth')), vDep='depth'; elseif any(strcmp(vnames, 'depth_std')), vDep='depth_std'; elseif any(strcmp(vnames, 'level')), vDep='level'; elseif any(strcmp(vnames, 'lev')), vDep='lev'; else, vDep='z'; end

Lon = ncread(SampleFile, vLon); Lat = ncread(SampleFile, vLat);
Depth = ncread(SampleFile, vDep); 
Lon = double(Lon); Lat = double(Lat);
if any(isnan(Lon)), Lon(isnan(Lon)) = 0; end
if any(isnan(Lat)), Lat(isnan(Lat)) = 0; end
Lon(Lon > 180) = Lon(Lon > 180) - 360;
Lat(Lat < -89.9) = -89.9; Lat(Lat > 89.9) = 89.9;
Depth = double(Depth);
if any(Depth < 0), Depth(Depth < 0) = 0; end

DepthIdx = find(Depth <= MaxDepth); Nx=length(Lon); Ny=length(Lat); Nz=length(DepthIdx);
fprintf('>> Grid: Nx=%d, Ny=%d, Nz=%d, MaxDepth=%.1f\n', Nx, Ny, Nz, Depth(DepthIdx(end)));

Sum_T = zeros(Nx, Ny, Nz); Sum_S = zeros(Nx, Ny, Nz); 
Cnt = zeros(Nx, Ny, Nz);

for y = Years
    try
        PatT = sprintf('temp.%d.nc', y);
        PatS = sprintf('sal.%d.nc', y);
        
        FT = fullfile(DataDir, 'Temperature', PatT);
        FS = fullfile(DataDir, 'Salinity', PatS);
        
        if ~exist(FT, 'file') || ~exist(FS, 'file')
            fprintf('   [Warning] File not found: %s or %s\n', PatT, PatS);
            continue;
        end

        vi = ncinfo(FT); vnames = {vi.Variables.Name};
        if any(strcmp(vnames, 'temp')), vT='temp'; elseif any(strcmp(vnames, 'temperature')), vT='temperature'; elseif any(strcmp(vnames, 'pt')), vT='pt'; else, vT='temp'; end
        
        vi = ncinfo(FS); vnames = {vi.Variables.Name};
        if any(strcmp(vnames, 'sal')), vS='sal'; elseif any(strcmp(vnames, 'salt')), vS='salt'; elseif any(strcmp(vnames, 'salinity')), vS='salinity'; elseif any(strcmp(vnames, 'so')), vS='so'; else, vS='sal'; end
        
        for m = 1:12
            T_raw = double(ncread(FT, vT, [1 1 1 m], [Inf Inf Nz 1]));
            S_raw = double(ncread(FS, vS, [1 1 1 m], [Inf Inf Nz 1]));
            
            mask = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
            
            T_raw(~mask) = 0;
            S_raw(~mask) = 0;
             
            Sum_T = Sum_T + T_raw;
            Sum_S = Sum_S + S_raw;
            Cnt = Cnt + double(mask);
        end
        fprintf('Year %d Done.\n', y);
        
    catch ME
        fprintf('   [Error] Year %d: %s\n', y, ME.message);
    end
end

Mean_T = Sum_T ./ Cnt;
Mean_S = Sum_S ./ Cnt;
Mean_T(Cnt==0) = NaN;
Mean_S(Cnt==0) = NaN;

fprintf('Mean_T: min=%.2f, max=%.2f, valid=%d, NaN=%d\n', ...
    min(Mean_T(:), [], 'omitnan'), max(Mean_T(:), [], 'omitnan'), sum(Cnt(:)>0), sum(Cnt(:)==0));
fprintf('Mean_S: min=%.2f, max=%.2f, valid=%d, NaN=%d\n', ...
    min(Mean_S(:), [], 'omitnan'), max(Mean_S(:), [], 'omitnan'), sum(Cnt(:)>0), sum(Cnt(:)==0));

SaveName = fullfile(OutputDir, 'Ishii_Mean_State.mat');
save(SaveName, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');
fprintf('Saved %s\n', SaveName);

warning(ws);
