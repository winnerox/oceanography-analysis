%% Calc_Mean_IAP.m
% =========================================================================
% 功能：计算 IAP 的平均态
% 输出：D:\work\IAP_TSLA_Terms\IAP_Mean_State.mat
% =========================================================================
clear; clc;

ws = warning('off', 'all');

DataDir = 'D:\work\IAP_05_24';
OutputDir = 'D:\work\IAP_TSLA_Terms';
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; MaxDepth = 2000;

List = dir(fullfile(DataDir, 'TEMP', '*.nc'));
SampleFile = fullfile(List(1).folder, List(1).name);

Lon = ncread(SampleFile, 'lon'); Lat = ncread(SampleFile, 'lat');
vi = ncinfo(SampleFile);
vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'depth_std')), vDepth='depth_std';
elseif any(strcmp(vnames, 'depth')), vDepth='depth'; 
elseif any(strcmp(vnames, 'level')), vDepth='level';
elseif any(strcmp(vnames, 'lev')), vDepth='lev';
else, error('Cannot find depth variable in %s', SampleFile); end

Depth = ncread(SampleFile, vDepth); 
DepthIdx = find(Depth <= MaxDepth); 
Nx = length(Lon); Ny = length(Lat); Nz = length(DepthIdx);

fprintf('Grid: %d x %d, Depth levels: %d (using %d <= %dm)\n', Nx, Ny, length(Depth), Nz, MaxDepth);

Sum_T = zeros(Nx, Ny, Nz); Sum_S = zeros(Nx, Ny, Nz); 
Cnt = zeros(Nx, Ny, Nz);

for y = Years
    for m = 1:12
        Pat = sprintf('*year_%d_month_%02d.nc', y, m);
        try
            sT = dir(fullfile(DataDir, 'TEMP', Pat));
            sS = dir(fullfile(DataDir, 'SALT', Pat));
            if isempty(sT) || isempty(sS), continue; end
            
            FT = fullfile(sT(1).folder, sT(1).name);
            FS = fullfile(sS(1).folder, sS(1).name);
            
            T_raw = ncread(FT, 'temp', [1 1 1], [Nz Inf Inf]);
            S_raw = ncread(FS, 'salinity', [1 1 1], [Nz Inf Inf]);
            
            T_raw = permute(double(T_raw), [2, 3, 1]);
            S_raw = permute(double(S_raw), [2, 3, 1]);
            
            mask = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
            
            T_raw(~mask) = 0;
            S_raw(~mask) = 0;
            
            Sum_T = Sum_T + T_raw;
            Sum_S = Sum_S + S_raw;
            Cnt = Cnt + double(mask);
        catch ME
            fprintf('Error reading %d-%02d: %s\n', y, m, ME.message);
        end
    end
    fprintf('Year %d Done.\n', y);
end

Mean_T = Sum_T ./ Cnt;
Mean_S = Sum_S ./ Cnt;
Mean_T(Cnt == 0) = NaN;
Mean_S(Cnt == 0) = NaN;

fprintf('Mean_T: min=%.2f, max=%.2f, valid=%d, NaN=%d\n', ...
    min(Mean_T(:), [], 'omitnan'), max(Mean_T(:), [], 'omitnan'), sum(Cnt(:)>0), sum(Cnt(:)==0));
fprintf('Mean_S: min=%.2f, max=%.2f, valid=%d, NaN=%d\n', ...
    min(Mean_S(:), [], 'omitnan'), max(Mean_S(:), [], 'omitnan'), sum(Cnt(:)>0), sum(Cnt(:)==0));

SaveName = fullfile(OutputDir, 'IAP_Mean_State.mat');
save(SaveName, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx', '-v7.3');
fprintf('Saved %s\n', SaveName);

warning(ws);
