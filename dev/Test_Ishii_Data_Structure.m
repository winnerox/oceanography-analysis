%% Test_Ishii_Data_Structure.m
% 检查 Ishii 数据的结构
clear; clc;

TFile = 'D:\work\Ishii_05_24\Temperature\temp.2005.nc';
SFile = 'D:\work\Ishii_05_24\Salinity\sal.2005.nc';

fprintf('========================================\n');
fprintf('  检查 Ishii 数据结构\n');
fprintf('========================================\n\n');

%% 1. 检查温度文件
fprintf('>>> 温度文件: %s\n', TFile);
warning('off', 'all');
InfoT = ncinfo(TFile);
warning('on', 'all');

fprintf('\n=== 变量列表 ===\n');
for i = 1:length(InfoT.Variables)
    v = InfoT.Variables(i);
    dims = '';
    for d = 1:length(v.Dimensions)
        dims = [dims, sprintf('%s(%d) ', v.Dimensions(d).Name, v.Dimensions(d).Length)];
    end
    fprintf('  %-15s: %s\n', v.Name, dims);
end

fprintf('\n=== 维度列表 ===\n');
for i = 1:length(InfoT.Dimensions)
    fprintf('  %-15s: %d\n', InfoT.Dimensions(i).Name, InfoT.Dimensions(i).Length);
end

%% 2. 尝试读取数据
fprintf('\n=== 读取数据测试 ===\n');

% 尝试常见的变量名
VarNamesT = {InfoT.Variables.Name};
TVar = '';
for v = {'temp', 'temperature', 'theta', 't', 'TEMP', 'Temperature'}
    if any(strcmpi(v{1}, VarNamesT))
        TVar = v{1};
        break;
    end
end

if isempty(TVar)
    fprintf('  温度变量: 未找到标准名，变量列表: %s\n', strjoin(VarNamesT, ', '));
else
    fprintf('  温度变量: %s\n', TVar);
    T = ncread(TFile, TVar);
    fprintf('  T 维度: %s\n', mat2str(size(T)));
    fprintf('  T 数据类型: %s\n', class(T));
    T = double(T);
    fprintf('  T 原始范围: [%.2f, %.2f]\n', min(T(:)), max(T(:)));
    % 清理 (常见 FillValue: 1e30, 999, -999)
    T(abs(T) > 1e10) = NaN;
    T(T > 100 | T < -10) = NaN;
    fprintf('  T 有效范围: [%.2f, %.2f]\n', min(T(:)), max(T(:)));
    fprintf('  T NaN 比例: %.1f%%\n', sum(isnan(T(:)))/numel(T)*100);
end

%% 3. 检查坐标变量
fprintf('\n=== 坐标变量 ===\n');
for v = {'lon', 'longitude', 'x', 'LON'}
    if any(strcmpi(v{1}, VarNamesT))
        lon = ncread(TFile, v{1});
        fprintf('  Lon (%s): [%.2f, %.2f], size=%d\n', v{1}, min(lon), max(lon), length(lon));
        break;
    end
end

for v = {'lat', 'latitude', 'y', 'LAT'}
    if any(strcmpi(v{1}, VarNamesT))
        lat = ncread(TFile, v{1});
        fprintf('  Lat (%s): [%.2f, %.2f], size=%d\n', v{1}, min(lat), max(lat), length(lat));
        break;
    end
end

for v = {'depth', 'lev', 'level', 'z', 'DEPTH', 'LEV'}
    if any(strcmpi(v{1}, VarNamesT))
        depth = ncread(TFile, v{1});
        fprintf('  Depth (%s): [%.1f, %.1f] m, size=%d\n', v{1}, min(depth), max(depth), length(depth));
        break;
    end
end

for v = {'time', 't', 'TIME', 'month'}
    if any(strcmpi(v{1}, VarNamesT))
        time = ncread(TFile, v{1});
        fprintf('  Time (%s): size=%d, values=%s\n', v{1}, length(time), mat2str(time(:)'));
        break;
    end
end

%% 4. 检查盐度文件
fprintf('\n>>> 盐度文件: %s\n', SFile);
warning('off', 'all');
InfoS = ncinfo(SFile);
warning('on', 'all');

VarNamesS = {InfoS.Variables.Name};
fprintf('  盐度变量列表: %s\n', strjoin(VarNamesS, ', '));

SVar = '';
for v = {'sal', 'salinity', 's', 'SALT', 'SAL', 'Salinity'}
    if any(strcmpi(v{1}, VarNamesS))
        SVar = v{1};
        break;
    end
end

if ~isempty(SVar)
    fprintf('  盐度变量: %s\n', SVar);
    S = ncread(SFile, SVar);
    fprintf('  S 维度: %s\n', mat2str(size(S)));
    S = double(S);
    fprintf('  S 原始范围: [%.2f, %.2f]\n', min(S(:)), max(S(:)));
    S(abs(S) > 1e10) = NaN;
    S(S > 100 | S < 0) = NaN;
    fprintf('  S 有效范围: [%.2f, %.2f]\n', min(S(:)), max(S(:)));
    fprintf('  S NaN 比例: %.1f%%\n', sum(isnan(S(:)))/numel(S)*100);
end

fprintf('\n>>> 测试完成! <<<\n');
