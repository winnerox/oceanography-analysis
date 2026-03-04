%% debug_en4_mean.m
YearDir = 'D:\work\EN4_analyses_c13_last20years\EN.4.2.2.analyses.c13.2005';
List = dir(fullfile(YearDir, '*.nc'));
f = fullfile(List(1).folder, List(1).name);
fprintf('File: %s\n', f);

Depth = ncread(f, 'depth');
fprintf('Depth: %d levels, max=%.1f\n', length(Depth), max(Depth));
DepthIdx = find(Depth <= 2000);
Nz = length(DepthIdx);
fprintf('Using %d levels <= 2000m\n', Nz);

T_raw = double(ncread(f, 'temperature', [1 1 1 1], [Inf Inf Nz 1]));
S_raw = double(ncread(f, 'salinity', [1 1 1 1], [Inf Inf Nz 1]));
fprintf('T_raw size: %s\n', mat2str(size(T_raw)));
fprintf('T_raw range before: %.2f to %.2f\n', min(T_raw(:), [], 'omitnan'), max(T_raw(:), [], 'omitnan'));

if mean(T_raw(:), 'omitnan') > 100
    T_raw = T_raw - 273.15;
    fprintf('Converted from Kelvin\n');
end
fprintf('T_raw range after: %.2f to %.2f\n', min(T_raw(:), [], 'omitnan'), max(T_raw(:), [], 'omitnan'));
fprintf('S_raw range: %.2f to %.2f\n', min(S_raw(:), [], 'omitnan'), max(S_raw(:), [], 'omitnan'));

mask = (T_raw > -100 & T_raw < 100) & (S_raw >= 0 & S_raw < 100);
fprintf('Valid points: %d / %d\n', sum(mask(:)), numel(mask));
