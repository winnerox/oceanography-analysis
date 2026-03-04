%% Calc_Mean_EN4.m
% =========================================================================
% 功能：计算 EN4 的平均态
% 输出：D:\work\EN4_TSLA_Terms\EN4_Mean_State.mat
% =========================================================================
clear; clc;

ws = warning('off', 'all');

DataDir = 'D:\work\EN4_analyses_c13_last20years';
OutputDir = 'D:\work\EN4_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; MaxDepth = 2000;

YearDir = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', Years(1)));
List = dir(fullfile(YearDir, '*.nc'));
if isempty(List), error('Cannot find EN4 data in %s', YearDir); end
SampleFile = fullfile(List(1).folder, List(1).name);

Lon = double(ncread(SampleFile, 'lon')); 
Lat = double(ncread(SampleFile, 'lat'));
vi = ncinfo(SampleFile);
vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth')), vDepth='depth'; 
elseif any(strcmp(vnames, 'deptht')), vDepth='deptht'; 
elseif any(strcmp(vnames, 'level')), vDepth='level';
elseif any(strcmp(vnames, 'lev')), vDepth='lev';
else, error('Cannot find depth variable in %s', SampleFile); end

Depth = double(ncread(SampleFile, vDepth)); 
DepthIdx = find(Depth <= MaxDepth); 
Nx = length(Lon); Ny = length(Lat); Nz = length(DepthIdx);

fprintf('Grid: %d x %d, Depth levels: %d\n', Nx, Ny, Nz);

Sum_T = zeros(Nx, Ny, Nz); Sum_S = zeros(Nx, Ny, Nz); 
Cnt = zeros(Nx, Ny, Nz);

for y = Years
    YearDir = fullfile(DataDir, sprintf('EN.4.2.2.analyses.c13.%d', y));
    if ~exist(YearDir, 'dir')
        fprintf('Year %d: Directory not found\n', y);
        continue;
    end
    
    for m = 1:12
        FilePat = sprintf('EN.4.2.2.f.analysis.c13.%d%02d.nc', y, m);
        try
            sFiles = dir(fullfile(YearDir, FilePat));
            if isempty(sFiles), continue; end
            F = fullfile(sFiles(1).folder, sFiles(1).name);
            
            vi = ncinfo(F); vnames = {vi.Variables.Name};
            if any(strcmp(vnames, 'temperature')), vT='temperature'; 
            elseif any(strcmp(vnames, 'temp')), vT='temp'; 
            elseif any(strcmp(vnames, 'thetao')), vT='thetao'; 
            else, vT=vnames{1}; end
            
            if any(strcmp(vnames, 'salinity')), vS='salinity'; 
            elseif any(strcmp(vnames, 'salt')), vS='salt'; 
            elseif any(strcmp(vnames, 'so')), vS='so'; 
            else, vS=vnames{2}; end
            
            T_raw = double(ncread(F, vT, [1 1 1 1], [Inf Inf Nz 1]));
            S_raw = double(ncread(F, vS, [1 1 1 1], [Inf Inf Nz 1]));
            
            if mean(T_raw(:), 'omitnan') > 100
                T_raw = T_raw - 273.15;
            end
            
            mask = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
            T_raw(~mask) = 0;
            S_raw(~mask) = 0;
            
            Sum_T = Sum_T + T_raw;
            Sum_S = Sum_S + S_raw;
            Cnt = Cnt + double(mask);
        catch
        end
    end
    fprintf('Year %d Done.\n', y);
end

Mean_T = Sum_T ./ Cnt;
Mean_S = Sum_S ./ Cnt;
Mean_T(Cnt == 0) = NaN;
Mean_S(Cnt == 0) = NaN;

fprintf('Mean_T: min=%.2f, max=%.2f\n', min(Mean_T(:), [], 'omitnan'), max(Mean_T(:), [], 'omitnan'));
fprintf('Mean_S: min=%.2f, max=%.2f\n', min(Mean_S(:), [], 'omitnan'), max(Mean_S(:), [], 'omitnan'));

SaveName = fullfile(OutputDir, 'EN4_Mean_State.mat');
save(SaveName, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx', '-v7.3');
fprintf('Saved %s\n', SaveName);

warning(ws);
