%% Explain_Cross_Terms.m
% =========================================================================
% 功能：解释混合项中包含的具体偏导项及其物理意义
% =========================================================================
clear; clc;

% 定义输出文件
output_file = 'Cross_Terms_Explanation.txt';
fid = fopen(output_file, 'w');

fprintf(fid, '=== 混合项 (Cross Terms) 偏导项详解 ===\n\n');
fprintf('=== 混合项 (Cross Terms) 偏导项详解 ===\n\n');

% 定义阶数范围
max_order = 8;

fprintf(fid, '混合项包含的偏导项 (2-%d阶):\n', max_order);
fprintf('混合项包含的偏导项 (2-%d阶):\n', max_order);
fprintf(fid, '=====================================\n\n');
fprintf('=====================================\n\n');

% 遍历每个阶数
for n = 2:max_order
    fprintf(fid, '第 %d 阶混合项:\n', n);
    fprintf('第 %d 阶混合项:\n', n);
    fprintf(fid, '-------------------------------------\n');
    fprintf('-------------------------------------\n');
    
    % 遍历所有可能的混合偏导数
    for k = 1:n-1
        n_T = n - k;
        n_S = k;
        
        % 计算二项式系数
        C = 1/(factorial(k) * factorial(n-k));
        
        % 生成偏导数描述
        if n_T == 1 && n_S == 1
            desc = '热膨胀-盐收缩交叉项';
        elseif n_T == 2 && n_S == 1
            desc = 'Cabbeling-盐收缩交叉项';
        elseif n_T == 1 && n_S == 2
            desc = '热膨胀-盐收缩非线性交叉项';
        else
            desc = sprintf('∂^%dρ/∂T^%d∂S^%d 混合偏导数', n, n_T, n_S);
        end
        
        fprintf(fid, '  - T^%d S^%d: %s (系数: %.6f)\n', n_T, n_S, desc, C);
        fprintf('  - T^%d S^%d: %s (系数: %.6f)\n', n_T, n_S, desc, C);
    end
    
    fprintf(fid, '\n');
    fprintf('\n');
end

fprintf(fid, '=== 物理意义说明 ===\n');
fprintf('=== 物理意义说明 ===\n');
fprintf(fid, '=====================================\n\n');
fprintf('=====================================\n\n');
fprintf(fid, '1. 混合项的物理意义：\n');
fprintf('1. 混合项的物理意义：\n');
fprintf(fid, '   - 捕捉温度和盐度变化对密度的非线性交叉影响\n');
fprintf('   - 捕捉温度和盐度变化对密度的非线性交叉影响\n');
fprintf(fid, '   - 描述温度和盐度异常如何共同影响海平面变化\n');
fprintf('   - 描述温度和盐度异常如何共同影响海平面变化\n');
fprintf(fid, '   - 高阶混合项反映了更复杂的热力学相互作用\n\n');
fprintf('   - 高阶混合项反映了更复杂的热力学相互作用\n\n');
fprintf(fid, '2. 关键偏导项解释：\n');
fprintf('2. 关键偏导项解释：\n');
fprintf(fid, '   - T^1 S^1: 一阶交叉项，反映温度和盐度变化的线性相互作用\n');
fprintf('   - T^1 S^1: 一阶交叉项，反映温度和盐度变化的线性相互作用\n');
fprintf(fid, '   - T^2 S^1: 温度二阶、盐度一阶交叉项，反映温度非线性与盐度的相互作用\n');
fprintf('   - T^2 S^1: 温度二阶、盐度一阶交叉项，反映温度非线性与盐度的相互作用\n');
fprintf(fid, '   - T^1 S^2: 温度一阶、盐度二阶交叉项，反映盐度非线性与温度的相互作用\n');
fprintf('   - T^1 S^2: 温度一阶、盐度二阶交叉项，反映盐度非线性与温度的相互作用\n');
fprintf(fid, '   - 高阶项: 反映更复杂的非线性热力学过程\n\n');
fprintf('   - 高阶项: 反映更复杂的非线性热力学过程\n\n');
fprintf(fid, '3. 计算方法：\n');
fprintf('3. 计算方法：\n');
fprintf(fid, '   - 对密度ρ关于温度(CT)和盐度(SA)进行混合偏导\n');
fprintf('   - 对密度ρ关于温度(CT)和盐度(SA)进行混合偏导\n');
fprintf(fid, '   - 使用二项式系数加权不同阶数的交叉项\n');
fprintf('   - 使用二项式系数加权不同阶数的交叉项\n');
fprintf(fid, '   - 乘以温度和盐度异常的相应幂次\n');
fprintf('   - 乘以温度和盐度异常的相应幂次\n');
fprintf(fid, '   - 对所有深度积分得到海平面贡献\n\n');
fprintf('   - 对所有深度积分得到海平面贡献\n\n');

fprintf(fid, '=== 数据文件结构 ===\n');
fprintf('=== 数据文件结构 ===\n');
fprintf(fid, '=====================================\n\n');
fprintf('=====================================\n\n');
fprintf(fid, '混合项文件中包含的变量：\n');
fprintf('混合项文件中包含的变量：\n');
fprintf(fid, '   - Cross_AllOrders: 混合项数据 [经度×纬度×时间×阶数]\n');
fprintf('   - Cross_AllOrders: 混合项数据 [经度×纬度×时间×阶数]\n');
fprintf(fid, '   - lon: 经度坐标\n');
fprintf('   - lon: 经度坐标\n');
fprintf(fid, '   - lat: 纬度坐标\n');
fprintf('   - lat: 纬度坐标\n');
fprintf(fid, '   - time_vec: 时间坐标\n\n');
fprintf('   - time_vec: 时间坐标\n\n');
fprintf(fid, '数据维度说明：\n');
fprintf('数据维度说明：\n');
fprintf(fid, '   - EN4: [360×173×240×8] (360经度×173纬度×240时间×8阶)\n');
fprintf('   - EN4: [360×173×240×8] (360经度×173纬度×240时间×8阶)\n');
fprintf(fid, '   - IAP: [360×180×240×8] (360经度×180纬度×240时间×8阶)\n\n');
fprintf('   - IAP: [360×180×240×8] (360经度×180纬度×240时间×8阶)\n\n');

fprintf(fid, '>>> 混合项偏导项解释完成! <<<\n');
fprintf('>>> 混合项偏导项解释完成! <<<\n');

fclose(fid);
fprintf('结果已保存到: %s\n', output_file);

