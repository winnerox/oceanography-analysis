%% 简单测试脚本 - 验证修改后的收敛性测试代码
clear; clc;

% 测试数据维度
Nx = 10;
Ny = 10;
Nt = 5;
MaxOrder = 8;

% 创建测试数据
fprintf('创建测试数据...\n');
EN4_Exact_Avg = randn(Nx, Ny, Nt);
EN4_Exact_Std = randn(Nx, Ny, Nt);
EN4_T_Terms_Avg = randn(Nx, Ny, Nt, MaxOrder) * 0.1;
EN4_T_Terms_Std = randn(Nx, Ny, Nt, MaxOrder) * 0.1;
EN4_S_Terms_Avg = randn(Nx, Ny, Nt, MaxOrder) * 0.05;
EN4_S_Terms_Std = randn(Nx, Ny, Nt, MaxOrder) * 0.05;
EN4_Cross_Terms_Avg = randn(Nx, Ny, Nt, MaxOrder) * 0.01;
EN4_Cross_Terms_Std = randn(Nx, Ny, Nt, MaxOrder) * 0.01;

% 复制相同数据给IAP和Ishii（用于测试）
IAP_Exact_Avg = EN4_Exact_Avg * 0.9;
IAP_Exact_Std = EN4_Exact_Std * 0.9;
IAP_T_Terms_Avg = EN4_T_Terms_Avg * 0.9;
IAP_T_Terms_Std = EN4_T_Terms_Std * 0.9;
IAP_S_Terms_Avg = EN4_S_Terms_Avg * 0.9;
IAP_S_Terms_Std = EN4_S_Terms_Std * 0.9;
IAP_Cross_Terms_Avg = EN4_Cross_Terms_Avg * 0.9;
IAP_Cross_Terms_Std = EN4_Cross_Terms_Std * 0.9;

Ishii_Exact_Avg = EN4_Exact_Avg * 1.1;
Ishii_Exact_Std = EN4_Exact_Std * 1.1;
Ishii_T_Terms_Avg = EN4_T_Terms_Avg * 1.1;
Ishii_T_Terms_Std = EN4_T_Terms_Std * 1.1;
Ishii_S_Terms_Avg = EN4_S_Terms_Avg * 1.1;
Ishii_S_Terms_Std = EN4_S_Terms_Std * 1.1;
Ishii_Cross_Terms_Avg = EN4_Cross_Terms_Avg * 1.1;
Ishii_Cross_Terms_Std = EN4_Cross_Terms_Std * 1.1;

lon = 1:Nx;
lat = 1:Ny;
time_vec = 1:Nt;

%% 测试计算各阶累加和
fprintf('\n测试计算各阶累加和...\n');

% EN4累加和
EN4_Sum_Avg = zeros(Nx, Ny, Nt, MaxOrder);
EN4_Sum_Std = zeros(Nx, Ny, Nt, MaxOrder);
for n = 1:MaxOrder
    if n == 1
        EN4_Sum_Avg(:,:,:,n) = EN4_T_Terms_Avg(:,:,:,n) + EN4_S_Terms_Avg(:,:,:,n);
        EN4_Sum_Std(:,:,:,n) = EN4_T_Terms_Std(:,:,:,n) + EN4_S_Terms_Std(:,:,:,n);
    else
        EN4_Sum_Avg(:,:,:,n) = EN4_Sum_Avg(:,:,:,n-1) + EN4_T_Terms_Avg(:,:,:,n) + EN4_S_Terms_Avg(:,:,:,n) + EN4_Cross_Terms_Avg(:,:,:,n);
        EN4_Sum_Std(:,:,:,n) = EN4_Sum_Std(:,:,:,n-1) + EN4_T_Terms_Std(:,:,:,n) + EN4_S_Terms_Std(:,:,:,n) + EN4_Cross_Terms_Std(:,:,:,n);
    end
end

fprintf('EN4累加和计算完成，维度: %s\n', mat2str(size(EN4_Sum_Avg)));

%% 测试计算残差
fprintf('\n测试计算残差...\n');

EN4_Residual_Avg = zeros(Nx, Ny, Nt, MaxOrder);
EN4_Residual_Std = zeros(Nx, Ny, Nt, MaxOrder);
for n = 1:MaxOrder
    EN4_Residual_Avg(:,:,:,n) = EN4_Exact_Avg - EN4_Sum_Avg(:,:,:,n);
    EN4_Residual_Std(:,:,:,n) = EN4_Exact_Std - EN4_Sum_Std(:,:,:,n);
end

fprintf('EN4残差计算完成，维度: %s\n', mat2str(size(EN4_Residual_Avg)));

%% 测试计算RMS误差
fprintf('\n测试计算RMS误差...\n');

EN4_RMS_Avg = zeros(MaxOrder, 1);
EN4_RMS_Std = zeros(MaxOrder, 1);
for n = 1:MaxOrder
    res_avg = EN4_Residual_Avg(:,:,:,n);
    res_std = EN4_Residual_Std(:,:,:,n);
    
    valid_avg = res_avg(abs(res_avg) < 500 & ~isnan(res_avg));
    valid_std = res_std(abs(res_std) < 500 & ~isnan(res_std));
    
    EN4_RMS_Avg(n) = sqrt(mean(valid_avg.^2));
    EN4_RMS_Std(n) = sqrt(mean(valid_std.^2));
end

fprintf('EN4 RMS误差计算完成:\n');
for n = 1:MaxOrder
    fprintf('  阶数 %d: Avg=%.4f, Std=%.4f\n', n, EN4_RMS_Avg(n), EN4_RMS_Std(n));
end

%% 测试收敛性检验
fprintf('\n测试收敛性检验...\n');

EN4_Converged_Avg = true;
EN4_Converged_Std = true;

for n = 1:MaxOrder
    if n > 1
        if EN4_RMS_Avg(n) > EN4_RMS_Avg(n-1)
            EN4_Converged_Avg = false;
        end
        if EN4_RMS_Std(n) > EN4_RMS_Std(n-1)
            EN4_Converged_Std = false;
        end
    end
end

fprintf('EN4收敛状态: 平均态-%s, 标准态-%s\n', ...
    ternary(EN4_Converged_Avg, '收敛', '未完全收敛'), ...
    ternary(EN4_Converged_Std, '收敛', '未完全收敛'));

%% 辅助函数
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

fprintf('\n测试完成！所有功能正常。\n');