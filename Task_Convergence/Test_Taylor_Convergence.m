%% Test_Taylor_Convergence.m
% =========================================================================
% 功能：测试泰勒展开收敛性 - 验证精确解与各阶累加和的差是否逐渐减少
% 扩展：增加对IAP和Ishii数据的收敛性测试
% 修复：修复 OOM 内存溢出、NaN 相加污染、精确解单位不一致(m vs mm)的问题
% =========================================================================
clear; clc;
addpath('D:\work');

% 定义数据目录
EN4_DataDir = 'D:\work\EN4_TSLA_Terms';
IAP_DataDir = 'D:\work\IAP_mat_data';
Ishii_DataDir = 'D:\work\Ishii_mat_data';

%% 1. 加载数据 (加入智能变量提取与单位对齐)
fprintf('加载数据...\n');

% 定义内部加载辅助函数
load_exact = @(filepath) get_exact_data(filepath);

% ========== EN4 数据 ==========
fprintf('加载 EN4 数据...\n');
% 精确解
EN4_Exact_Avg = load_exact(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Avg.mat'));
EN4_Exact_Std = load_exact(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Std.mat'));

% 提取坐标
tmp_coord = load(fullfile(EN4_DataDir, 'EN4_Formula11_Exact_Avg.mat'), 'lon', 'lat', 'time_vec');
lon = tmp_coord.lon; lat = tmp_coord.lat; time_vec = tmp_coord.time_vec;

% 展开项
load(fullfile(EN4_DataDir, 'EN4_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders'); EN4_T_Terms_Avg = TSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders'); EN4_T_Terms_Std = TSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders'); EN4_S_Terms_Avg = SSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders'); EN4_S_Terms_Std = SSLA_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders'); EN4_Cross_Terms_Avg = Cross_AllOrders;
load(fullfile(EN4_DataDir, 'EN4_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders'); EN4_Cross_Terms_Std = Cross_AllOrders;

% ========== IAP 数据 ==========
fprintf('加载 IAP 数据...\n');
% 精确解
IAP_Exact_Avg = load_exact(fullfile(IAP_DataDir, 'IAP_Formula11_Exact_Avg.mat'));
IAP_Exact_Std = load_exact(fullfile(IAP_DataDir, 'IAP_Formula11_Exact_StdRef.mat'));

% 展开项
load(fullfile(IAP_DataDir, 'IAP_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders'); IAP_T_Terms_Avg = TSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders'); IAP_T_Terms_Std = TSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders'); IAP_S_Terms_Avg = SSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders'); IAP_S_Terms_Std = SSLA_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders'); IAP_Cross_Terms_Avg = Cross_AllOrders;
load(fullfile(IAP_DataDir, 'IAP_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders'); IAP_Cross_Terms_Std = Cross_AllOrders;

% ========== Ishii 数据 ==========
fprintf('加载 Ishii 数据...\n');
% 精确解
Ishii_Exact_Avg = load_exact(fullfile(Ishii_DataDir, 'Ishii_Formula11_Exact_Avg.mat'));
Ishii_Exact_Std = load_exact(fullfile(Ishii_DataDir, 'Ishii_Formula11_Exact_StdRef.mat'));

% 展开项
load(fullfile(Ishii_DataDir, 'Ishii_TSLA_Terms_1to8_Average.mat'), 'TSLA_AllOrders'); Ishii_T_Terms_Avg = TSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_TSLA_Terms_1to8_StdRef.mat'), 'TSLA_AllOrders'); Ishii_T_Terms_Std = TSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_S_Terms_1to8_Average.mat'), 'SSLA_AllOrders'); Ishii_S_Terms_Avg = SSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_S_Terms_1to8_StdRef.mat'), 'SSLA_AllOrders'); Ishii_S_Terms_Std = SSLA_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_Cross_Terms_1to8_Average.mat'), 'Cross_AllOrders'); Ishii_Cross_Terms_Avg = Cross_AllOrders;
load(fullfile(Ishii_DataDir, 'Ishii_Cross_Terms_1to8_StdRef.mat'), 'Cross_AllOrders'); Ishii_Cross_Terms_Std = Cross_AllOrders;

% 获取各数据集维度
[EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder] = size(EN4_T_Terms_Avg);
[IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder] = size(IAP_T_Terms_Avg);
[Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder] = size(Ishii_T_Terms_Avg);

fprintf('EN4 数据维度: %d x %d x %d, 阶数: %d\n', EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder);
fprintf('IAP 数据维度: %d x %d x %d, 阶数: %d\n', IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder);
fprintf('Ishii 数据维度: %d x %d x %d, 阶数: %d\n', Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder);

% 使用EN4的阶数作为最大阶数（假设所有数据集阶数相同）
MaxOrder = EN4_MaxOrder;

%% 2. 计算各阶累加和 (使用安全防 NaN 累加并使用 single 节约内存)
fprintf('计算各阶累加和...\n');

% ========== EN4 累加和 ==========
fprintf('计算 EN4 累加和...\n');
EN4_Sum_Avg = zeros(EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder, 'single');
EN4_Sum_Std = zeros(EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder, 'single');
for n = 1:EN4_MaxOrder
    if n == 1
        EN4_Sum_Avg(:,:,:,n) = safe_nansum(EN4_T_Terms_Avg(:,:,:,n), EN4_S_Terms_Avg(:,:,:,n));
        EN4_Sum_Std(:,:,:,n) = safe_nansum(EN4_T_Terms_Std(:,:,:,n), EN4_S_Terms_Std(:,:,:,n));
    else
        EN4_Sum_Avg(:,:,:,n) = safe_nansum(EN4_Sum_Avg(:,:,:,n-1), EN4_T_Terms_Avg(:,:,:,n), EN4_S_Terms_Avg(:,:,:,n), EN4_Cross_Terms_Avg(:,:,:,n));
        EN4_Sum_Std(:,:,:,n) = safe_nansum(EN4_Sum_Std(:,:,:,n-1), EN4_T_Terms_Std(:,:,:,n), EN4_S_Terms_Std(:,:,:,n), EN4_Cross_Terms_Std(:,:,:,n));
    end
end

% ========== IAP 累加和 ==========
fprintf('计算 IAP 累加和...\n');
IAP_Sum_Avg = zeros(IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder, 'single');
IAP_Sum_Std = zeros(IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder, 'single');
for n = 1:IAP_MaxOrder
    if n == 1
        IAP_Sum_Avg(:,:,:,n) = safe_nansum(IAP_T_Terms_Avg(:,:,:,n), IAP_S_Terms_Avg(:,:,:,n));
        IAP_Sum_Std(:,:,:,n) = safe_nansum(IAP_T_Terms_Std(:,:,:,n), IAP_S_Terms_Std(:,:,:,n));
    else
        IAP_Sum_Avg(:,:,:,n) = safe_nansum(IAP_Sum_Avg(:,:,:,n-1), IAP_T_Terms_Avg(:,:,:,n), IAP_S_Terms_Avg(:,:,:,n), IAP_Cross_Terms_Avg(:,:,:,n));
        IAP_Sum_Std(:,:,:,n) = safe_nansum(IAP_Sum_Std(:,:,:,n-1), IAP_T_Terms_Std(:,:,:,n), IAP_S_Terms_Std(:,:,:,n), IAP_Cross_Terms_Std(:,:,:,n));
    end
end

% ========== Ishii 累加和 ==========
fprintf('计算 Ishii 累加和...\n');
Ishii_Sum_Avg = zeros(Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder, 'single');
Ishii_Sum_Std = zeros(Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder, 'single');
for n = 1:Ishii_MaxOrder
    if n == 1
        Ishii_Sum_Avg(:,:,:,n) = safe_nansum(Ishii_T_Terms_Avg(:,:,:,n), Ishii_S_Terms_Avg(:,:,:,n));
        Ishii_Sum_Std(:,:,:,n) = safe_nansum(Ishii_T_Terms_Std(:,:,:,n), Ishii_S_Terms_Std(:,:,:,n));
    else
        Ishii_Sum_Avg(:,:,:,n) = safe_nansum(Ishii_Sum_Avg(:,:,:,n-1), Ishii_T_Terms_Avg(:,:,:,n), Ishii_S_Terms_Avg(:,:,:,n), Ishii_Cross_Terms_Avg(:,:,:,n));
        Ishii_Sum_Std(:,:,:,n) = safe_nansum(Ishii_Sum_Std(:,:,:,n-1), Ishii_T_Terms_Std(:,:,:,n), Ishii_S_Terms_Std(:,:,:,n), Ishii_Cross_Terms_Std(:,:,:,n));
    end
end

%% 3. 计算残差
fprintf('计算残差...\n');

% ========== EN4 残差 ==========
fprintf('计算 EN4 残差...\n');
EN4_Residual_Avg = zeros(EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder, 'single');
EN4_Residual_Std = zeros(EN4_Nx, EN4_Ny, EN4_Nt, EN4_MaxOrder, 'single');
for n = 1:EN4_MaxOrder
    EN4_Residual_Avg(:,:,:,n) = EN4_Exact_Avg - EN4_Sum_Avg(:,:,:,n);
    EN4_Residual_Std(:,:,:,n) = EN4_Exact_Std - EN4_Sum_Std(:,:,:,n);
end

% ========== IAP 残差 ==========
fprintf('计算 IAP 残差...\n');
IAP_Residual_Avg = zeros(IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder, 'single');
IAP_Residual_Std = zeros(IAP_Nx, IAP_Ny, IAP_Nt, IAP_MaxOrder, 'single');
for n = 1:IAP_MaxOrder
    IAP_Residual_Avg(:,:,:,n) = IAP_Exact_Avg - IAP_Sum_Avg(:,:,:,n);
    IAP_Residual_Std(:,:,:,n) = IAP_Exact_Std - IAP_Sum_Std(:,:,:,n);
end

% ========== Ishii 残差 ==========
fprintf('计算 Ishii 残差...\n');
Ishii_Residual_Avg = zeros(Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder, 'single');
Ishii_Residual_Std = zeros(Ishii_Nx, Ishii_Ny, Ishii_Nt, Ishii_MaxOrder, 'single');
for n = 1:Ishii_MaxOrder
    Ishii_Residual_Avg(:,:,:,n) = Ishii_Exact_Avg - Ishii_Sum_Avg(:,:,:,n);
    Ishii_Residual_Std(:,:,:,n) = Ishii_Exact_Std - Ishii_Sum_Std(:,:,:,n);
end

%% 4. 计算全局RMS误差 (加入异常点过滤)
fprintf('计算全局RMS误差...\n');

% ========== EN4 RMS误差 ==========
fprintf('计算 EN4 RMS误差...\n');
EN4_RMS_Avg = zeros(EN4_MaxOrder, 1);
EN4_RMS_Std = zeros(EN4_MaxOrder, 1);
for n = 1:EN4_MaxOrder
    res_avg = EN4_Residual_Avg(:,:,:,n);
    res_std = EN4_Residual_Std(:,:,:,n);
    
    valid_avg = res_avg(abs(res_avg) < 500 & ~isnan(res_avg));
    valid_std = res_std(abs(res_std) < 500 & ~isnan(res_std));
    
    EN4_RMS_Avg(n) = sqrt(mean(valid_avg.^2));
    EN4_RMS_Std(n) = sqrt(mean(valid_std.^2));
end

% ========== IAP RMS误差 ==========
fprintf('计算 IAP RMS误差...\n');
IAP_RMS_Avg = zeros(IAP_MaxOrder, 1);
IAP_RMS_Std = zeros(IAP_MaxOrder, 1);
for n = 1:IAP_MaxOrder
    res_avg = IAP_Residual_Avg(:,:,:,n);
    res_std = IAP_Residual_Std(:,:,:,n);
    
    valid_avg = res_avg(abs(res_avg) < 500 & ~isnan(res_avg));
    valid_std = res_std(abs(res_std) < 500 & ~isnan(res_std));
    
    IAP_RMS_Avg(n) = sqrt(mean(valid_avg.^2));
    IAP_RMS_Std(n) = sqrt(mean(valid_std.^2));
end

% ========== Ishii RMS误差 ==========
fprintf('计算 Ishii RMS误差...\n');
Ishii_RMS_Avg = zeros(Ishii_MaxOrder, 1);
Ishii_RMS_Std = zeros(Ishii_MaxOrder, 1);
for n = 1:Ishii_MaxOrder
    res_avg = Ishii_Residual_Avg(:,:,:,n);
    res_std = Ishii_Residual_Std(:,:,:,n);
    
    valid_avg = res_avg(abs(res_avg) < 500 & ~isnan(res_avg));
    valid_std = res_std(abs(res_std) < 500 & ~isnan(res_std));
    
    Ishii_RMS_Avg(n) = sqrt(mean(valid_avg.^2));
    Ishii_RMS_Std(n) = sqrt(mean(valid_std.^2));
end

%% 5. 检验收敛性
fprintf('\n========== EN4 收敛性检验 ==========\n');
fprintf('阶数 | RMS_Avg (mm) | RMS_Std (mm) | 收敛比率\n');
fprintf('-----|--------------|--------------|--------\n');

EN4_Converged_Avg = true;
EN4_Converged_Std = true;

for n = 1:EN4_MaxOrder
    if n == 1
        fprintf('  %d  |   %8.4f   |   %8.4f   |   -\n', n, EN4_RMS_Avg(n), EN4_RMS_Std(n));
    else
        ratio_Avg = EN4_RMS_Avg(n-1) / EN4_RMS_Avg(n);
        ratio_Std = EN4_RMS_Std(n-1) / EN4_RMS_Std(n);
        fprintf('  %d  |   %8.4f   |   %8.4f   |  %.2f / %.2f\n', ...
            n, EN4_RMS_Avg(n), EN4_RMS_Std(n), ratio_Avg, ratio_Std);
        
        if EN4_RMS_Avg(n) > EN4_RMS_Avg(n-1)
            EN4_Converged_Avg = false;
        end
        if EN4_RMS_Std(n) > EN4_RMS_Std(n-1)
            EN4_Converged_Std = false;
        end
    end
end

fprintf('\n========== IAP 收敛性检验 ==========\n');
fprintf('阶数 | RMS_Avg (mm) | RMS_Std (mm) | 收敛比率\n');
fprintf('-----|--------------|--------------|--------\n');

IAP_Converged_Avg = true;
IAP_Converged_Std = true;

for n = 1:IAP_MaxOrder
    if n == 1
        fprintf('  %d  |   %8.4f   |   %8.4f   |   -\n', n, IAP_RMS_Avg(n), IAP_RMS_Std(n));
    else
        ratio_Avg = IAP_RMS_Avg(n-1) / IAP_RMS_Avg(n);
        ratio_Std = IAP_RMS_Std(n-1) / IAP_RMS_Std(n);
        fprintf('  %d  |   %8.4f   |   %8.4f   |  %.2f / %.2f\n', ...
            n, IAP_RMS_Avg(n), IAP_RMS_Std(n), ratio_Avg, ratio_Std);
        
        if IAP_RMS_Avg(n) > IAP_RMS_Avg(n-1)
            IAP_Converged_Avg = false;
        end
        if IAP_RMS_Std(n) > IAP_RMS_Std(n-1)
            IAP_Converged_Std = false;
        end
    end
end

fprintf('\n========== Ishii 收敛性检验 ==========\n');
fprintf('阶数 | RMS_Avg (mm) | RMS_Std (mm) | 收敛比率\n');
fprintf('-----|--------------|--------------|--------\n');

Ishii_Converged_Avg = true;
Ishii_Converged_Std = true;

for n = 1:Ishii_MaxOrder
    if n == 1
        fprintf('  %d  |   %8.4f   |   %8.4f   |   -\n', n, Ishii_RMS_Avg(n), Ishii_RMS_Std(n));
    else
        ratio_Avg = Ishii_RMS_Avg(n-1) / Ishii_RMS_Avg(n);
        ratio_Std = Ishii_RMS_Std(n-1) / Ishii_RMS_Std(n);
        fprintf('  %d  |   %8.4f   |   %8.4f   |  %.2f / %.2f\n', ...
            n, Ishii_RMS_Avg(n), Ishii_RMS_Std(n), ratio_Avg, ratio_Std);
        
        if Ishii_RMS_Avg(n) > Ishii_RMS_Avg(n-1)
            Ishii_Converged_Avg = false;
        end
        if Ishii_RMS_Std(n) > Ishii_RMS_Std(n-1)
            Ishii_Converged_Std = false;
        end
    end
end

%% 6. 输出结论
fprintf('\n========== EN4 结论 ==========\n');
if EN4_Converged_Avg
    fprintf('✅ EN4 平均态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ EN4 平均态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

if EN4_Converged_Std
    fprintf('✅ EN4 标准态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ EN4 标准态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

fprintf('\n========== IAP 结论 ==========\n');
if IAP_Converged_Avg
    fprintf('✅ IAP 平均态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ IAP 平均态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

if IAP_Converged_Std
    fprintf('✅ IAP 标准态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ IAP 标准态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

fprintf('\n========== Ishii 结论 ==========\n');
if Ishii_Converged_Avg
    fprintf('✅ Ishii 平均态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ Ishii 平均态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

if Ishii_Converged_Std
    fprintf('✅ Ishii 标准态：泰勒展开收敛，残差随阶数增加而减小\n');
else
    fprintf('❌ Ishii 标准态：泰勒展开未完全收敛，存在阶数残差增大\n');
end

%% 7. 绘图
% 图1：三个数据集的平均态收敛性对比
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
semilogy(1:EN4_MaxOrder, EN4_RMS_Avg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(1:IAP_MaxOrder, IAP_RMS_Avg, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(1:Ishii_MaxOrder, Ishii_RMS_Avg, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('平均态收敛性对比', 'FontSize', 14);
legend('EN4', 'IAP', 'Ishii', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

subplot(2,2,2);
semilogy(1:EN4_MaxOrder, EN4_RMS_Std, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(1:IAP_MaxOrder, IAP_RMS_Std, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(1:Ishii_MaxOrder, Ishii_RMS_Std, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('标准态收敛性对比', 'FontSize', 14);
legend('EN4', 'IAP', 'Ishii', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

subplot(2,2,3);
bar_data = [EN4_RMS_Avg, IAP_RMS_Avg, Ishii_RMS_Avg];
bar(bar_data);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('平均态各阶残差对比', 'FontSize', 14);
legend('EN4', 'IAP', 'Ishii', 'Location', 'best');
set(gca, 'FontSize', 11);

subplot(2,2,4);
bar_data = [EN4_RMS_Std, IAP_RMS_Std, Ishii_RMS_Std];
bar(bar_data);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('标准态各阶残差对比', 'FontSize', 14);
legend('EN4', 'IAP', 'Ishii', 'Location', 'best');
set(gca, 'FontSize', 11);

sgtitle('三套数据泰勒展开收敛性检验', 'FontSize', 16);

% 图2：各数据集独立绘图
figure('Position', [100, 100, 1500, 400]);

subplot(1,3,1);
semilogy(1:EN4_MaxOrder, EN4_RMS_Avg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(1:EN4_MaxOrder, EN4_RMS_Std, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('EN4 收敛性', 'FontSize', 14);
legend('平均态', '标准态', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

subplot(1,3,2);
semilogy(1:IAP_MaxOrder, IAP_RMS_Avg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(1:IAP_MaxOrder, IAP_RMS_Std, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('IAP 收敛性', 'FontSize', 14);
legend('平均态', '标准态', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

subplot(1,3,3);
semilogy(1:Ishii_MaxOrder, Ishii_RMS_Avg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(1:Ishii_MaxOrder, Ishii_RMS_Std, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('泰勒展开阶数', 'FontSize', 12);
ylabel('RMS残差 (mm)', 'FontSize', 12);
title('Ishii 收敛性', 'FontSize', 14);
legend('平均态', '标准态', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

sgtitle('各数据集独立收敛性检验', 'FontSize', 16);

%% 8. 保存结果
% 保存EN4结果
EN4_ResultFile = fullfile(EN4_DataDir, 'Taylor_Convergence_Test.mat');
save(EN4_ResultFile, 'EN4_RMS_Avg', 'EN4_RMS_Std', 'EN4_Residual_Avg', 'EN4_Residual_Std', ...
    'EN4_Sum_Avg', 'EN4_Sum_Std', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('\nEN4结果已保存至: %s\n', EN4_ResultFile);

% 保存IAP结果
IAP_ResultFile = fullfile(IAP_DataDir, 'Taylor_Convergence_Test.mat');
save(IAP_ResultFile, 'IAP_RMS_Avg', 'IAP_RMS_Std', 'IAP_Residual_Avg', 'IAP_Residual_Std', ...
    'IAP_Sum_Avg', 'IAP_Sum_Std', '-v7.3');
fprintf('IAP结果已保存至: %s\n', IAP_ResultFile);

% 保存Ishii结果
Ishii_ResultFile = fullfile(Ishii_DataDir, 'Taylor_Convergence_Test.mat');
save(Ishii_ResultFile, 'Ishii_RMS_Avg', 'Ishii_RMS_Std', 'Ishii_Residual_Avg', 'Ishii_Residual_Std', ...
    'Ishii_Sum_Avg', 'Ishii_Sum_Std', '-v7.3');
fprintf('Ishii结果已保存至: %s\n', Ishii_ResultFile);

% 保存综合结果
Combined_ResultFile = fullfile('D:\work\Task_Convergence', 'Combined_Taylor_Convergence_Test.mat');
save(Combined_ResultFile, 'EN4_RMS_Avg', 'EN4_RMS_Std', 'IAP_RMS_Avg', 'IAP_RMS_Std', ...
    'Ishii_RMS_Avg', 'Ishii_RMS_Std', 'lon', 'lat', 'time_vec', '-v7.3');
fprintf('综合结果已保存至: %s\n', Combined_ResultFile);

%% 9. 详细统计
fprintf('\n========== 详细统计 ==========\n');

% EN4统计
fprintf('\n--- EN4 统计 ---\n');
fprintf('平均态最终残差: %.4f mm (%d阶)\n', EN4_RMS_Avg(EN4_MaxOrder), EN4_MaxOrder);
fprintf('标准态最终残差: %.4f mm (%d阶)\n', EN4_RMS_Std(EN4_MaxOrder), EN4_MaxOrder);
fprintf('平均态收敛比率 (%d阶/1阶): %.2f%%\n', EN4_MaxOrder, (EN4_RMS_Avg(EN4_MaxOrder)/EN4_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (%d阶/1阶): %.2f%%\n', EN4_MaxOrder, (EN4_RMS_Std(EN4_MaxOrder)/EN4_RMS_Std(1))*100);

% IAP统计
fprintf('\n--- IAP 统计 ---\n');
fprintf('平均态最终残差: %.4f mm (%d阶)\n', IAP_RMS_Avg(IAP_MaxOrder), IAP_MaxOrder);
fprintf('标准态最终残差: %.4f mm (%d阶)\n', IAP_RMS_Std(IAP_MaxOrder), IAP_MaxOrder);
fprintf('平均态收敛比率 (%d阶/1阶): %.2f%%\n', IAP_MaxOrder, (IAP_RMS_Avg(IAP_MaxOrder)/IAP_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (%d阶/1阶): %.2f%%\n', IAP_MaxOrder, (IAP_RMS_Std(IAP_MaxOrder)/IAP_RMS_Std(1))*100);

% Ishii统计
fprintf('\n--- Ishii 统计 ---\n');
fprintf('平均态最终残差: %.4f mm (%d阶)\n', Ishii_RMS_Avg(Ishii_MaxOrder), Ishii_MaxOrder);
fprintf('标准态最终残差: %.4f mm (%d阶)\n', Ishii_RMS_Std(Ishii_MaxOrder), Ishii_MaxOrder);
fprintf('平均态收敛比率 (%d阶/1阶): %.2f%%\n', Ishii_MaxOrder, (Ishii_RMS_Avg(Ishii_MaxOrder)/Ishii_RMS_Avg(1))*100);
fprintf('标准态收敛比率 (%d阶/1阶): %.2f%%\n', Ishii_MaxOrder, (Ishii_RMS_Std(Ishii_MaxOrder)/Ishii_RMS_Std(1))*100);

% 综合对比
fprintf('\n--- 三套数据对比 ---\n');
fprintf('平均态最终残差对比 (EN4/IAP/Ishii): %.4f / %.4f / %.4f mm\n', ...
    EN4_RMS_Avg(EN4_MaxOrder), IAP_RMS_Avg(IAP_MaxOrder), Ishii_RMS_Avg(Ishii_MaxOrder));
fprintf('标准态最终残差对比 (EN4/IAP/Ishii): %.4f / %.4f / %.4f mm\n', ...
    EN4_RMS_Std(EN4_MaxOrder), IAP_RMS_Std(IAP_MaxOrder), Ishii_RMS_Std(Ishii_MaxOrder));

% 收敛性总结
fprintf('\n--- 收敛性总结 ---\n');
fprintf('EN4 收敛状态: 平均态-%s, 标准态-%s\n', ...
    ternary(EN4_Converged_Avg, '收敛', '未完全收敛'), ...
    ternary(EN4_Converged_Std, '收敛', '未完全收敛'));
fprintf('IAP 收敛状态: 平均态-%s, 标准态-%s\n', ...
    ternary(IAP_Converged_Avg, '收敛', '未完全收敛'), ...
    ternary(IAP_Converged_Std, '收敛', '未完全收敛'));
fprintf('Ishii 收敛状态: 平均态-%s, 标准态-%s\n', ...
    ternary(Ishii_Converged_Avg, '收敛', '未完全收敛'), ...
    ternary(Ishii_Converged_Std, '收敛', '未完全收敛'));

%% 10. 辅助函数
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function data = get_exact_data(file_path)
    % 智能加载 Exact 数据，规避变量名不一致问题，并检测/对齐单位到毫米(mm)
    tmp = load(file_path);
    fns = fieldnames(tmp);
    for i = 1:length(fns)
        if contains(fns{i}, 'TSLA') && isnumeric(tmp.(fns{i}))
            data = single(tmp.(fns{i}));
            
            % 智能判断单位：如果绝对值的均值非常小 (<0.5)，大概率单位是米
            if mean(abs(data(:)), 'omitnan') < 0.5
                fprintf('   [单位修复] 检测到 %s 可能是 m，已自动 * 1000 转换为 mm\n', fns{i});
                data = data * 1000;
            end
            return;
        end
    end
    error('未在 %s 中找到包含 TSLA 的有效数据矩阵', file_path);
end

function res = safe_nansum(varargin)
    % 针对多维数组的安全相加函数 (继承 single 类型以节约内存)
    res = zeros(size(varargin{1}), class(varargin{1}));
    all_nan = true(size(varargin{1})); 
    
    for k = 1:nargin
        mat = varargin{k};
        all_nan = all_nan & isnan(mat);
        mat(isnan(mat)) = 0; % 仅在计算加法时视为 0
        res = res + mat;
    end
    
    % 如果某处在所有项中都是 NaN(比如陆地)，则恢复其 NaN 状态
    res(all_nan) = NaN;
end