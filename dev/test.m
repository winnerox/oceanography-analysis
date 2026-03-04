% 快速检查 IAP 原始数据质量
TFile = 'D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc';
SFile = 'D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc';

% 读取
warning('off', 'all');
T = ncread(TFile, 'temp');
S = ncread(SFile, 'salinity');
warning('on', 'all');

% 检查
fprintf('T 维度: %s\n', mat2str(size(T)));
fprintf('S 维度: %s\n', mat2str(size(S)));
fprintf('T 数据类型: %s\n', class(T));
fprintf('S 数据类型: %s\n', class(S));

% 原始范围 (包含 FillValue)
fprintf('T 原始范围: [%.2f, %.2f]\n', min(T(:)), max(T(:)));
fprintf('S 原始范围: [%.2f, %.2f]\n', min(S(:)), max(S(:)));

% 清理后范围
T(T > 100 | T < -10) = NaN;
S(S > 100 | S < 0) = NaN;
fprintf('T 有效范围: [%.2f, %.2f]\n', min(T(:)), max(T(:)));
fprintf('S 有效范围: [%.2f, %.2f]\n', min(S(:)), max(S(:)));

% NaN 比例
fprintf('T NaN 比例: %.1f%%\n', sum(isnan(T(:)))/numel(T)*100);
fprintf('S NaN 比例: %.1f%%\n', sum(isnan(S(:)))/numel(S)*100);