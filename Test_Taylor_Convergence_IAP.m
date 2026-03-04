%% Check_Residuals_Convergence.m
% =========================================================================
% 功能：计算并可视化海平面高阶泰勒展开的“残差收敛情况”
% 分析项：精确解 vs (1阶), (1+2阶), (1+2+3阶) ... (1~8阶累加)
% =========================================================================
clear; clc; close all;

%% 1. 配置参数 (根据需要修改此处)
DataDir  = 'D:\work\EN4_mat_data'; % 数据所在文件夹
Dataset  = 'EN4';                  % 数据集名称: 'EN4', 'IAP', 或 'Ishii'
Mode     = 'StdRef';               % 模式: 'StdRef' (标准态) 或 'Average' (平均态)

%% 2. 自动构建文件名并加载数据
fprintf('>> [1/4] 正在加载 %s (%s) 的精确解与泰勒项数据...\n', Dataset, Mode);

if strcmp(Mode, 'StdRef')
    File_Exact = fullfile(DataDir, sprintf('%s_Formula11_Exact_StdRef.mat', Dataset));
    File_T     = fullfile(DataDir, sprintf('%s_TSLA_Terms_1to8_StdRef.mat', Dataset));
    File_S     = fullfile(DataDir, sprintf('%s_S_Terms_1to8_StdRef.mat', Dataset));
    File_C     = fullfile(DataDir, sprintf('%s_Cross_Terms_1to8_StdRef.mat', Dataset));
else
    File_Exact = fullfile(DataDir, sprintf('%s_Formula11_Exact_Avg.mat', Dataset));
    File_T     = fullfile(DataDir, sprintf('%s_TSLA_Terms_1to8_Average.mat', Dataset));
    File_S     = fullfile(DataDir, sprintf('%s_S_Terms_1to8_Average.mat', Dataset));
    File_C     = fullfile(DataDir, sprintf('%s_Cross_Terms_1to8_Average.mat', Dataset));
end

try
    % 按模式精准加载对应的变量，彻底消灭未找到变量的警告
    if strcmp(Mode, 'StdRef')
        load(File_Exact, 'TSLA_Exact_StdRef', 'lat', 'lon', 'Time_Axis');
        Exact_3D = TSLA_Exact_StdRef;
    else
        load(File_Exact, 'TSLA_Exact_Avg', 'lat', 'lon', 'Time_Axis');
        Exact_3D = TSLA_Exact_Avg;
    end
    
    load(File_T, 'TSLA_AllOrders');
    load(File_S, 'SSLA_AllOrders');
    load(File_C, 'Cross_AllOrders');
catch ME
    error('❌ 数据加载失败，请检查路径和文件名是否正确: %s', ME.message);
end

[Nx, Ny, Nt, MaxOrder] = size(TSLA_AllOrders);
fprintf('   ✅ 加载成功！时间步: %d, 最大阶数: %d\n', Nt, MaxOrder);

%% 3. 计算累加的泰勒展开总和 (Cumulative Sum)
fprintf('>> [2/4] 计算各阶累加总和与 3D 残差...\n');

Taylor_CumSum_3D = zeros(Nx, Ny, Nt, MaxOrder, 'single');
Residuals_3D     = zeros(Nx, Ny, Nt, MaxOrder, 'single');

% 1阶单独计算 (由于交叉项从2阶开始，1阶的C项全为NaN或0)
% 将 NaN 视为 0 进行安全的加法
T1 = TSLA_AllOrders(:,:,:,1); T1(isnan(T1)) = 0;
S1 = SSLA_AllOrders(:,:,:,1); S1(isnan(S1)) = 0;
Taylor_CumSum_3D(:,:,:,1) = T1 + S1;

% 生成用于恢复真实 NaN 陆地掩膜的基准
LandMask = isnan(Exact_3D);

for k = 2:MaxOrder
    Tk = TSLA_AllOrders(:,:,:,k); Tk(isnan(Tk)) = 0;
    Sk = SSLA_AllOrders(:,:,:,k); Sk(isnan(Sk)) = 0;
    Ck = Cross_AllOrders(:,:,:,k); Ck(isnan(Ck)) = 0;
    
    % 当前 k 阶累加 = 前 k-1 阶累加 + 第 k 阶的三项
    Taylor_CumSum_3D(:,:,:,k) = Taylor_CumSum_3D(:,:,:,k-1) + Tk + Sk + Ck;
end

% 把陆地掩膜重新贴回去，并计算 3D 残差
for k = 1:MaxOrder
    tmp_sum = Taylor_CumSum_3D(:,:,:,k);
    tmp_sum(LandMask) = NaN;
    Taylor_CumSum_3D(:,:,:,k) = tmp_sum;
    
    Residuals_3D(:,:,:,k) = Exact_3D - tmp_sum;
end

%% 4. 面积加权计算全球平均 (GMSL)
fprintf('>> [3/4] 面积加权计算全球平均时序 (转换单位为 mm)...\n');

[LAT, ~] = meshgrid(lat, lon);
Weight = cosd(LAT); % 面积加权系数

GMSL_Exact    = zeros(Nt, 1);
GMSL_Taylor   = zeros(Nt, MaxOrder);
GMSL_Residual = zeros(Nt, MaxOrder);

for t = 1:Nt
    % 获取当前月的有效水柱掩膜
    exact_slice = Exact_3D(:,:,t);
    valid = ~isnan(exact_slice);
    
    W = Weight(valid);
    Sum_W = sum(W);
    
    % 精确解的全球平均 (转为 mm)
    GMSL_Exact(t) = sum(exact_slice(valid) .* W) / Sum_W * 1000;
    
    for k = 1:MaxOrder
        taylor_slice = Taylor_CumSum_3D(:,:,t,k);
        res_slice    = Residuals_3D(:,:,t,k);
        
        GMSL_Taylor(t, k)   = sum(taylor_slice(valid) .* W) / Sum_W * 1000;
        GMSL_Residual(t, k) = sum(res_slice(valid) .* W) / Sum_W * 1000;
    end
end

% 计算各阶残差的时序均方根误差 (RMS) 和 标准差 (Std)
RMS_Residual = sqrt(mean(GMSL_Residual.^2, 1));
STD_Residual = std(GMSL_Residual, 0, 1);

%% 4.5 在控制台输出收敛性检验表
fprintf('\n========== %s 收敛性检验 (%s) ==========\n', Dataset, Mode);
fprintf('阶数 | RMS_Avg (mm) | RMS_Std (mm) | 收敛比率\n');
fprintf('-----|--------------|--------------|----------------\n');

Converged_Avg = true;
Converged_Std = true;

for n = 1:MaxOrder
    if n == 1
        fprintf('  %d  |   %8.4f   |   %8.4f   |   -\n', n, RMS_Residual(n), STD_Residual(n));
    else
        ratio_Avg = RMS_Residual(n-1) / RMS_Residual(n);
        ratio_Std = STD_Residual(n-1) / STD_Residual(n);
        fprintf('  %d  |   %8.4f   |   %8.4f   | %7.2f / %5.2f\n', ...
            n, RMS_Residual(n), STD_Residual(n), ratio_Avg, ratio_Std);
        
        if RMS_Residual(n) > RMS_Residual(n-1)
            Converged_Avg = false;
        end
        if STD_Residual(n) > STD_Residual(n-1)
            Converged_Std = false;
        end
    end
end

fprintf('----------------------------------------------------\n');
if Converged_Avg && Converged_Std
    fprintf('>>> 结论: %s 模型在 %d 阶内呈现严格的单调收敛！\n', Dataset, MaxOrder);
else
    fprintf('>>> ⚠️ 注意: 模型在高阶可能存在截断误差发散或浮点极限波动。\n');
end
fprintf('====================================================\n\n');

%% 5. 极客级可视化
fprintf('>> [4/4] 绘制收敛情况图表...\n');

% -----------------------------------------------------------
% 图 1：海平面时序对比 (Exact vs Order 1 vs Order 2 vs Order 3)
% -----------------------------------------------------------
figure('Name', 'GMSL Time Series', 'Position', [100, 100, 900, 400]);
plot(Time_Axis, GMSL_Exact, 'k', 'LineWidth', 2, 'DisplayName', 'Exact Solution'); hold on;
plot(Time_Axis, GMSL_Taylor(:,1), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Order 1 (Linear)');
plot(Time_Axis, GMSL_Taylor(:,2), 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Order \leq 2');
plot(Time_Axis, GMSL_Taylor(:,3), 'g:', 'LineWidth', 1.5, 'DisplayName', 'Order \leq 3');
grid on; box on;
xlabel('Year'); ylabel('Sea Level Anomaly (mm)');
title(sprintf('[%s - %s] Global Mean Sea Level Anomaly', Dataset, Mode));
legend('Location', 'best');
set(gca, 'FontSize', 12);

% -----------------------------------------------------------
% 图 2：残差时序衰减 (Residuals over time)
% -----------------------------------------------------------
figure('Name', 'Residuals Time Series', 'Position', [150, 150, 900, 400]);
plot(Time_Axis, GMSL_Residual(:,1), 'b', 'LineWidth', 1.5, 'DisplayName', 'Residual (1st Order)'); hold on;
plot(Time_Axis, GMSL_Residual(:,2), 'r', 'LineWidth', 1.5, 'DisplayName', 'Residual (up to 2nd)');
plot(Time_Axis, GMSL_Residual(:,3), 'g', 'LineWidth', 1.5, 'DisplayName', 'Residual (up to 3rd)');
plot(Time_Axis, GMSL_Residual(:,4), 'm', 'LineWidth', 1.5, 'DisplayName', 'Residual (up to 4th)');
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
grid on; box on;
xlabel('Year'); ylabel('Residual Error (mm)');
title(sprintf('[%s - %s] Residual Decay (Exact - Taylor Sum)', Dataset, Mode));
legend('Location', 'best');
set(gca, 'FontSize', 12);

% -----------------------------------------------------------
% 图 3：对数收敛瀑布图 (最核心的结论图)
% -----------------------------------------------------------
figure('Name', 'Convergence Waterfall', 'Position', [200, 200, 600, 500]);
semilogy(1:MaxOrder, RMS_Residual, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on; box on;
set(gca, 'YScale', 'log'); % 强制对数坐标
xticks(1:MaxOrder);
xlabel('Maximum Taylor Order (k)');
ylabel('RMS Residual Error (mm)');
title(sprintf('[%s] Exponential Convergence of Taylor Series', Dataset));
% 添加数值标签
for i = 1:MaxOrder
    text(i, RMS_Residual(i), sprintf('  %.1e', RMS_Residual(i)), 'VerticalAlignment', 'bottom', 'FontSize', 10);
end
set(gca, 'FontSize', 12);

fprintf('   🎉 残差分析完成！请查看弹出的 3 张分析图表。\n');