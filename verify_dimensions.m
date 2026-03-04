%% 验证维度一致性修复
clear; clc;
addpath('D:\work');

%% 1. 测试数据读取
fprintf('测试数据读取和维度转换...\n');
DataRoot = 'D:\work\IAP_05_24';
TempDir = fullfile(DataRoot, 'TEMP');
SaltDir = fullfile(DataRoot, 'SALT');

FileT = fullfile(TempDir, 'IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc');
FileS = fullfile(SaltDir, 'IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc');

% 读取坐标
Lon = single(ncread(FileT, 'lon'));
Lat = single(ncread(FileT, 'lat'));
Depth = single(ncread(FileT, 'depth_std'));

% 统一 lon 坐标
Lon = mod(Lon, 360);

fprintf('坐标维度:\n');
fprintf('  Lon: %d 点 (%.1f 到 %.1f)\n', length(Lon), min(Lon), max(Lon));
fprintf('  Lat: %d 点 (%.1f 到 %.1f)\n', length(Lat), min(Lat), max(Lat));
fprintf('  Depth: %d 层\n', length(Depth));

%% 2. 测试数据读取和 permute
MaxDepth = 2000;
DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);

% 读取原始数据
temp_raw = ncread(FileT, 'temp');
salt_raw = ncread(FileS, 'salinity');

fprintf('\n原始数据维度:\n');
fprintf('  temp_raw: %s [depth, lon, lat]\n', mat2str(size(temp_raw)));
fprintf('  salt_raw: %s [depth, lon, lat]\n', mat2str(size(salt_raw)));

% 应用 permute
temp_permuted = permute(temp_raw, [3, 1, 2]);  % [lat, depth, lon]
salt_permuted = permute(salt_raw, [3, 1, 2]);  % [lat, depth, lon]

fprintf('\npermute 后数据维度:\n');
fprintf('  temp_permuted: %s [lat, depth, lon]\n', mat2str(size(temp_permuted)));
fprintf('  salt_permuted: %s [lat, depth, lon]\n', mat2str(size(salt_permuted)));

%% 3. 生成网格并验证维度一致性
fprintf('\n生成网格...\n');
[LAT_3D, DEPTH_3D, LON_3D] = ndgrid(Lat, CalcDepth, Lon);
P_3D = single(gsw_p_from_z(-DEPTH_3D, LAT_3D));
[Ny, Nz, Nx] = size(P_3D);

fprintf('网格维度:\n');
fprintf('  LAT_3D: %s [lat, depth, lon]\n', mat2str(size(LAT_3D)));
fprintf('  DEPTH_3D: %s [lat, depth, lon]\n', mat2str(size(DEPTH_3D)));
fprintf('  LON_3D: %s [lat, depth, lon]\n', mat2str(size(LON_3D)));
fprintf('  P_3D: %s [lat, depth, lon]\n', mat2str(size(P_3D)));

% 截取 permuted 数据以匹配 CalcDepth
temp_cropped = temp_permuted(:, DepthIdx, :);
salt_cropped = salt_permuted(:, DepthIdx, :);

fprintf('\n截取后的数据维度 (匹配 CalcDepth):\n');
fprintf('  temp_cropped: %s [lat, depth, lon]\n', mat2str(size(temp_cropped)));
fprintf('  salt_cropped: %s [lat, depth, lon]\n', mat2str(size(salt_cropped)));

%% 4. 维度一致性检查
fprintf('\n维度一致性检查:\n');
dim_match_temp = isequal(size(temp_cropped), size(P_3D));
dim_match_salt = isequal(size(salt_cropped), size(P_3D));

fprintf('  temp_cropped 与 P_3D 维度匹配: %s\n', mat2str(dim_match_temp));
fprintf('  salt_cropped 与 P_3D 维度匹配: %s\n', mat2str(dim_match_salt));

if dim_match_temp && dim_match_salt
    fprintf('✅ 所有数据维度完全一致！\n');
else
    fprintf('❌ 维度不一致！需要进一步调试。\n');
    
    % 显示详细维度信息
    fprintf('\n详细维度对比:\n');
    fprintf('  temp_cropped: %s\n', mat2str(size(temp_cropped)));
    fprintf('  P_3D: %s\n', mat2str(size(P_3D)));
    fprintf('  差异: %s\n', mat2str(size(temp_cropped) - size(P_3D)));
end

%% 5. 测试索引访问
fprintf('\n测试索引访问...\n');
test_lat_idx = 90;   % 赤道附近
test_depth_idx = 10; % 第10层
test_lon_idx = 180;  % 0度经度

% 从网格获取坐标
test_lat = LAT_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_lon = LON_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_depth = DEPTH_3D(test_lat_idx, test_depth_idx, test_lon_idx);
test_p = P_3D(test_lat_idx, test_depth_idx, test_lon_idx);

% 从数据获取值
test_temp = temp_cropped(test_lat_idx, test_depth_idx, test_lon_idx);
test_salt = salt_cropped(test_lat_idx, test_depth_idx, test_lon_idx);

fprintf('  测试点 (lat_idx=%d, depth_idx=%d, lon_idx=%d):\n', ...
    test_lat_idx, test_depth_idx, test_lon_idx);
fprintf('    坐标: Lat=%.2f°, Lon=%.2f°, Depth=%.1f m, P=%.1f dbar\n', ...
    test_lat, test_lon, test_depth, test_p);
fprintf('    数据: Temp=%.2f°C, Salt=%.2f psu\n', test_temp, test_salt);

%% 6. 测试 GSW 函数
fprintf('\n测试 GSW 函数...\n');
try
    SA_test = gsw_SA_from_SP(double(test_salt), double(test_p), double(test_lon), double(test_lat));
    CT_test = gsw_CT_from_pt(SA_test, double(test_temp));
    
    fprintf('  GSW 转换成功:\n');
    fprintf('    绝对盐度: %.4f g/kg\n', SA_test);
    fprintf('    保守温度: %.4f°C\n', CT_test);
    
    % 测试密度计算
    rho_test = gsw_rho(SA_test, CT_test, test_p);
    fprintf('    海水密度: %.4f kg/m³\n', rho_test);
    
    fprintf('✅ GSW 函数测试通过！\n');
catch ME
    fprintf('❌ GSW 函数测试失败: %s\n', ME.message);
end

fprintf('\n✅ 维度验证完成！\n');