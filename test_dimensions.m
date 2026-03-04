%% 测试维度一致性
clear; clc;
addpath('D:\work');

%% 1. 读取样本文件
DataRoot = 'D:\work\IAP_05_24';
TempDir = fullfile(DataRoot, 'TEMP');
SampleFile = fullfile(TempDir, 'IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc');

if ~exist(SampleFile, 'file'), error('❌ 找不到样本文件'); end

Lon = single(ncread(SampleFile, 'lon'));
Lat = single(ncread(SampleFile, 'lat'));
Depth = single(ncread(SampleFile, 'depth_std'));

fprintf('原始网格维度:\n');
fprintf('  Lon: %d 点 (%.1f 到 %.1f)\n', length(Lon), min(Lon), max(Lon));
fprintf('  Lat: %d 点 (%.1f 到 %.1f)\n', length(Lat), min(Lat), max(Lat));
fprintf('  Depth: %d 层\n', length(Depth));

%% 2. 计算 CalcDepth
MaxDepth = 2000;
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
fprintf('\n计算深度 (<= %d m): %d 层\n', MaxDepth, length(CalcDepth));

%% 3. 生成网格
[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Ny, Nz, Nx] = size(P_3D);

fprintf('\n网格维度:\n');
fprintf('  LAT_3D: %d x %d x %d [lat, depth, lon]\n', size(LAT_3D));
fprintf('  DEPTH_3D: %d x %d x %d [lat, depth, lon]\n', size(DEPTH_3D));
fprintf('  LON_3D: %d x %d x %d [lat, depth, lon]\n', size(LON_3D));
fprintf('  P_3D: %d x %d x %d [lat, depth, lon]\n', size(P_3D));

%% 4. 测试读取数据
T_data = ncread(SampleFile, 'temp', [1 1 1], [length(Depth) Inf Inf]);
fprintf('\nNetCDF 读取维度:\n');
fprintf('  T_data (原始): %d x %d x %d [depth, lat, lon]\n', size(T_data));

T_permuted = permute(T_data, [2, 1, 3]);
fprintf('  T_permuted: %d x %d x %d [lat, depth, lon]\n', size(T_permuted));

%% 5. 验证维度一致性
fprintf('\n维度一致性检查:\n');
fprintf('  LAT_3D 与 T_permuted 维度匹配: %s\n', ...
    mat2str(isequal(size(LAT_3D), size(T_permuted))));
fprintf('  P_3D 与 T_permuted 维度匹配: %s\n', ...
    mat2str(isequal(size(P_3D), size(T_permuted))));

%% 6. 测试索引
fprintf('\n索引测试:\n');
test_lat_idx = 90;  % 赤道附近
test_lon_idx = 180; % 0度经度
test_depth_idx = 10; % 第10层

% 从网格中获取坐标
test_lat = LAT_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_lon = LON_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_depth = DEPTH_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_p = P_3D(test_lat_idx, test_depth_idx, test_lon_idx);

fprintf('  测试点 (lat_idx=%d, depth_idx=%d, lon_idx=%d):\n', ...
    test_lat_idx, test_depth_idx, test_lon_idx);
fprintf('    Lat: %.2f°, Lon: %.2f°, Depth: %.1f m, Pressure: %.1f dbar\n', ...
    test_lat, test_lon, test_depth, test_p);

%% 7. 测试 GSW 函数
fprintf('\nGSW 函数测试:\n');
test_temp = double(T_permuted(test_lat_idx, test_depth_idx, test_lon_idx));
test_sal = 35.0; % 假设盐度

SA_test = gsw_SA_from_SP(test_sal, test_p, test_lon, test_lat);
CT_test = gsw_CT_from_pt(SA_test, test_temp);

fprintf('  温度: %.2f°C, 盐度: %.2f psu\n', test_temp, test_sal);
fprintf('  绝对盐度: %.4f g/kg, 保守温度: %.4f°C\n', SA_test, CT_test);

fprintf('\n✅ 维度测试完成！\n');