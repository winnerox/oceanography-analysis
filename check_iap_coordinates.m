%% 检查IAP温度和盐度文件的经度坐标差异
clear; clc;

TempFile = 'D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc';
SaltFile = 'D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc';

fprintf('=== 检查IAP温度文件坐标 ===\n');
TempLon = ncread(TempFile, 'lon');
TempLat = ncread(TempFile, 'lat');
fprintf('温度经度: %d 点, 范围 %.2f 到 %.2f\n', length(TempLon), min(TempLon), max(TempLon));
fprintf('温度纬度: %d 点, 范围 %.2f 到 %.2f\n', length(TempLat), min(TempLat), max(TempLat));
fprintf('温度经度前5个值: %s\n', mat2str(TempLon(1:5)', 4));
fprintf('温度经度后5个值: %s\n', mat2str(TempLon(end-4:end)', 4));

fprintf('\n=== 检查IAP盐度文件坐标 ===\n');
SaltLon = ncread(SaltFile, 'lon');
SaltLat = ncread(SaltFile, 'lat');
fprintf('盐度经度: %d 点, 范围 %.2f 到 %.2f\n', length(SaltLon), min(SaltLon), max(SaltLon));
fprintf('盐度纬度: %d 点, 范围 %.2f 到 %.2f\n', length(SaltLat), min(SaltLat), max(SaltLat));
fprintf('盐度经度前5个值: %s\n', mat2str(SaltLon(1:5)', 4));
fprintf('盐度经度后5个值: %s\n', mat2str(SaltLon(end-4:end)', 4));

fprintf('\n=== 坐标差异分析 ===\n');
if length(TempLon) == length(SaltLon)
    lon_diff = TempLon - SaltLon;
    fprintf('经度差值 (前5个): %s\n', mat2str(lon_diff(1:5)', 4));
    fprintf('经度差值 (后5个): %s\n', mat2str(lon_diff(end-4:end)', 4));
    fprintf('最大经度偏差: %.4f 度\n', max(abs(lon_diff)));
else
    fprintf('警告: 经度点数不同! 温度=%d, 盐度=%d\n', length(TempLon), length(SaltLon));
end

if length(TempLat) == length(SaltLat)
    lat_diff = TempLat - SaltLat;
    fprintf('纬度差值 (前5个): %s\n', mat2str(lat_diff(1:5)', 4));
    fprintf('最大纬度偏差: %.4f 度\n', max(abs(lat_diff)));
else
    fprintf('警告: 纬度点数不同! 温度=%d, 盐度=%d\n', length(TempLat), length(SaltLat));
end

fprintf('\n=== 结论 ===\n');
if max(abs(lon_diff)) > 0.01
    fprintf('❌ 发现严重问题: 经度坐标存在 %.2f 度偏移!\n', max(abs(lon_diff)));
    fprintf('   必须对盐度数据进行空间插值对齐!\n');
else
    fprintf('✅ 经度坐标一致\n');
end
