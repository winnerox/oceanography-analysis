%% Get_Thermosteric_Poly_5th.m
% 功能: 生成 5 阶热比容密度近似函数 & 海平面变化公式
% 场景: 基础状态 SA=35, CT=10, p=1000
clear; clc; close all;

% 1. 启动引擎
engine = TEOS10_HighOrder_Engine();

% 2. 设定展开中心 (Base State)
SA_0 = 35.0;  CT_0 = 10.0;  p_0  = 1000;
H_depth = 2000;   % 假设水深 2000m
rho_0   = 1025.0; % 参考密度

fprintf('\n================================================================\n');
fprintf('   TEOS-10 热比容 5 阶近似公式生成器\n');
fprintf('   展开中心: SA=%.1f g/kg, CT=%.1f C\n', SA_0, CT_0);
fprintf('================================================================\n');

% 3. 计算 1-5 阶系数
% 泰勒公式: f(x) = c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5
% 其中 c_n = (1/n!) * f^(n)(0)
coeffs = zeros(5, 1);

fprintf('\n>>> [A] 密度变化多项式系数 (Delta_rho):\n');
fprintf('    Delta_rho = C1*dT + C2*dT^2 + C3*dT^3 + C4*dT^4 + C5*dT^5\n');
fprintf('----------------------------------------------------------------\n');

for n = 1:5
    % 计算 n 阶导数 (对 T 求 n 次, 对 S 求 0 次)
    deriv = engine.calculate_mixed(SA_0, CT_0, p_0, n, 0);
    
    % 计算系数
    c_n = deriv / factorial(n);
    coeffs(n) = c_n;
    
    % 打印系数
    if n==1, name='(线性膨胀)';
    elseif n==2, name='(非线性/Cabbeling)';
    else, name=sprintf('(%d阶修正)', n); end
    
    fprintf('    C%d = %15.6e  %s\n', n, c_n, name);
end

% 4. 转换为海平面变化公式系数
% Delta_h = - (H / rho_0) * Delta_rho
%         = K1*dT + K2*dT^2 + ...
k_coeffs = - (H_depth / rho_0) * coeffs * 1000; % 结果单位: mm

fprintf('\n>>> [B] 海平面上升多项式系数 (Delta_h, 单位: mm):\n');
fprintf('    Delta_h(mm) = K1*dT + K2*dT^2 + K3*dT^3 + K4*dT^4 + K5*dT^5\n');
fprintf('----------------------------------------------------------------\n');
for n = 1:5
    fprintf('    K%d = %15.6f mm/C^%d\n', n, k_coeffs(n), n);
end
fprintf('----------------------------------------------------------------\n');

% 5. 构造并验证函数
rho_poly = @(dT) coeffs(1)*dT + coeffs(2)*dT.^2 + coeffs(3)*dT.^3 + coeffs(4)*dT.^4 + coeffs(5)*dT.^5;
h_poly   = @(dT) k_coeffs(1)*dT + k_coeffs(2)*dT.^2 + k_coeffs(3)*dT.^3 + k_coeffs(4)*dT.^4 + k_coeffs(5)*dT.^5;

% 绘图验证 (模拟大温差 0 -> 15度)
dT_range = linspace(0, 15, 100);
h_vals = h_poly(dT_range);

% 计算各项贡献 (用于堆叠图)
term1 = k_coeffs(1)*dT_range;
term2 = k_coeffs(2)*dT_range.^2;
term3 = k_coeffs(3)*dT_range.^3;
term4 = k_coeffs(4)*dT_range.^4;
term5 = k_coeffs(5)*dT_range.^5;

figure('Color', 'w', 'Position', [100, 100, 900, 500]);

% 左图: 总预测曲线
subplot(1, 2, 1);
plot(dT_range, h_vals, 'r-', 'LineWidth', 2);
hold on;
plot(dT_range, term1, 'b--', 'LineWidth', 1.5); % 仅线性
grid on;
xlabel('温升 \Delta T (^\circC)');
ylabel('海平面上升 \Delta h (mm)');
title('5阶近似预测曲线 (2000m水深)');
legend('5阶全量预测', '仅1阶线性预测', 'Location', 'best');

% 右图: 高阶项的修正量 (2-5阶)
subplot(1, 2, 2);
plot(dT_range, term2, 'g-', 'LineWidth', 2, 'DisplayName', '2阶 (Cabbeling)');
hold on;
plot(dT_range, term3, 'm-', 'LineWidth', 2, 'DisplayName', '3阶修正');
plot(dT_range, term4, 'c-', 'LineWidth', 2, 'DisplayName', '4阶修正');
plot(dT_range, term5, 'k-', 'LineWidth', 2, 'DisplayName', '5阶修正');
grid on;
xlabel('温升 \Delta T (^\circC)');
ylabel('修正贡献量 (mm)');
title('各高阶项的修正贡献');
legend('Location', 'best');

% 6. 输出直接可用的公式文本
fprintf('\n>>> [C] 您的 5 阶近似公式 (可以直接复制到论文/代码):\n');
fprintf('Delta_rho(dT) = (%.3e)*dT + (%.3e)*dT^2 + (%.3e)*dT^3 + (%.3e)*dT^4 + (%.3e)*dT^5\n', ...
    coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5));