%% Debug_IAP_T3.m
clear; clc;

% 1. 加载 IAP 结果
DataFile = 'D:\work\MAT_Data\IAP_TSLA_Terms_1to8_Average.mat';
if ~exist(DataFile, 'file'), error('File not found'); end
fprintf('Loading %s...\n', DataFile);
load(DataFile, 'TSLA_AllOrders');

% 2. 统计各阶量级
fprintf('\n--- IAP Term Statistics (Mean Abs Value, mm) ---\n');
for i = 1:8
    data = TSLA_AllOrders(:,:,:,i);
    val = nanmean(abs(data(:)));
    fprintf('T%d: %.4f mm\n', i, val);
end

% 3. 检查 T1 vs T3 比率
T1 = TSLA_AllOrders(:,:,:,1);
T3 = TSLA_AllOrders(:,:,:,3);

Ratio = abs(T3) ./ abs(T1);
MeanRatio = nanmean(Ratio(:));
fprintf('\nMean Ratio |T3|/|T1|: %.4f (Expected < 0.05)\n', MeanRatio);

% 4. 尝试重算一个点的系数来验证 Engine
fprintf('\n--- Testing TEOS10 Engine for one point ---\n');
try
    Engine = TEOS10_General_Engine();
    % 假设 S=35, T=0, P=0 (表层标准态)
    SA = 35.16504; % gsw_SA_from_SP(35, 0, 0, 0)
    CT = 0;
    P = 0;
    
    Derivs = Engine.calculate_all_orders(SA, CT, P, 3);
    fprintf('Derivs (Sample):\n');
    disp(Derivs);
    
    % 手动估算 T3 贡献
    % Term3 = 1/6 * d3_T * dT^3
    % 假设 dT = 10度
    dT = 10;
    Term1 = Derivs.d1_T * dT;
    Term3 = (1/6) * Derivs.d3_T * (dT^3);
    
    fprintf('For dT = 10 deg:\n');
    fprintf('  Term1 ~ %.4e\n', Term1);
    fprintf('  Term3 ~ %.4e\n', Term3);
    fprintf('  Ratio T3/T1 ~ %.4f\n', abs(Term3/Term1));
    
catch ME
    fprintf('Engine Error: %s\n', ME.message);
end
