function gsw_rho_ttt_exact_test()
% GSW_RHO_TTT_EXACT_TEST
% 功能：
% 1. 自动从 GSW 源码中提取 75 项系数 (确保与官方版本一致)
% 2. 使用解析法 (链式法则) 计算密度对保守温度的三阶导数
% 3. 与高精度差分法进行对比验证
    
    clc; clear;
    
    %% 1. 准备测试数据
    % 建议先用转换函数处理原始数据，这里直接使用转换后的标准变量
    SA = 35.0;   % 绝对盐度 [g/kg]
    CT = 20.0;   % 保守温度 [deg C]
    p  = 1000.0; % 压力 [dbar]
    
    fprintf('======================================================\n');
    fprintf('   密度三阶导数计算：解析法 vs 差分法\n');
    fprintf('======================================================\n');
    fprintf('测试点: SA=%.2f, CT=%.2f, p=%.1f\n', SA, CT, p);

    %% 2. 获取系数 (核心步骤)
    % 尝试从 gsw_specvol.m 提取，如果失败则报错
    try
        [C, inds] = extract_gsw_coeffs();
        fprintf('>> [状态] 成功从 gsw_specvol.m 提取出 %d 个系数。\n', length(C));
    catch ME
        error('系数提取失败: %s\n请确保 gsw_specvol.m 在路径中且未被加密(.p)。', ME.message);
    end

    %% 3. 方法A: 解析求导 (Analytical)
    tic;
    [rho_ttt_exact, v_val, v_t, v_tt, v_ttt] = calc_rho_ttt_analytical(SA, CT, p, C, inds);
    time_exact = toc;

    %% 4. 方法B: 差分验证 (Finite Difference)
    % 使用 h=0.1 的中心差分作为"真值"
    tic;
    h = 0.1;
    r_m2 = gsw_rho(SA, CT-2*h, p);
    r_m1 = gsw_rho(SA, CT-h,   p);
    r_p1 = gsw_rho(SA, CT+h,   p);
    r_p2 = gsw_rho(SA, CT+2*h, p);
    rho_ttt_fd = (-0.5*r_m2 + 1.0*r_m1 - 1.0*r_p1 + 0.5*r_p2) / (h^3);
    time_fd = toc;

    %% 5. 结果对比
    fprintf('\n--- 最终结果 ---\n');
    fprintf('解析解 (Exact): %.12e (耗时 %.2f ms)\n', rho_ttt_exact, time_exact*1000);
    fprintf('差分解 (FD)   : %.12e (耗时 %.2f ms)\n', rho_ttt_fd, time_fd*1000);
    
    abs_err = abs(rho_ttt_exact - rho_ttt_fd);
    rel_err = abs_err / abs(rho_ttt_fd);
    
    fprintf('\n绝对误差      : %.2e\n', abs_err);
    fprintf('相对误差      : %.2e\n', rel_err);
    
    if rel_err < 1e-6
        fprintf('>> [通过] 验证成功！解析解与差分结果吻合。\n');
        
        % 额外输出比容的各阶导数供检查
        fprintf('\n--- 中间变量检查 (比容导数) ---\n');
        fprintf('v    : %.5e\n', v_val);
        fprintf('v_t  : %.5e\n', v_t);
        fprintf('v_tt : %.5e\n', v_tt);
        fprintf('v_ttt: %.5e\n', v_ttt);
    else
        fprintf('>> [警告] 存在差异，请检查系数提取是否完整。\n');
    end
    fprintf('======================================================\n');
end

%% ============================================================
%  核心函数 1: 链式解析求导
% ============================================================
function [rho_ttt, v, v_t, v_tt, v_ttt] = calc_rho_ttt_analytical(SA, CT, p, C, inds)
    % 1. 变量缩放 (Roquet et al., 2015)
    g_scl = 35.16504;
    t_scl = 40.0;
    p_scl = 10000.0;
    
    x = sqrt(max(0, SA) ./ g_scl); % 注意: 盐度是开根号的缩放变量
    y = CT ./ t_scl;
    z = p  ./ p_scl;
    
    % 2. 初始化累加变量
    v = 0;      % 0阶
    vy = 0;     % 对 y 的1阶
    vyy = 0;    % 对 y 的2阶
    vyyy = 0;   % 对 y 的3阶
    
    % 3. 遍历 75 项多项式求导
    for k = 1:length(C)
        coeff = C(k);
        i = inds(k, 1); % x指数 (SA)
        j = inds(k, 2); % y指数 (CT) - 关键求导项
        l = inds(k, 3); % z指数 (p)
        
        % 基础项 (与温度无关的部分)
        base = coeff * (x^i) * (z^l);
        
        % 0阶: v
        v = v + base * (y^j);
        
        % 1阶: dv/dy
        if j >= 1
            vy = vy + base * j * (y^(j-1));
        end
        
        % 2阶: d2v/dy2
        if j >= 2
            vyy = vyy + base * j*(j-1) * (y^(j-2));
        end
        
        % 3阶: d3v/dy3
        if j >= 3
            vyyy = vyyy + base * j*(j-1)*(j-2) * (y^(j-3));
        end
    end
    
    % 4. 恢复物理量导数 (Un-scaling)
    % 链式法则: d/dCT = (dy/dCT) * d/dy = (1/t_scl) * d/dy
    v_t   = vy   / t_scl;
    v_tt  = vyy  / (t_scl^2);
    v_ttt = vyyy / (t_scl^3);
    
    % 5. 从比容导数转为密度导数
    % rho = 1/v
    rho = 1/v;
    % 公式: rho''' = -v'''/v^2 + 6*v'*v''/v^3 - 6*(v')^3/v^4
    rho_ttt = -v_ttt * (rho^2) + 6 * v_t * v_tt * (rho^3) - 6 * (v_t^3) * (rho^4);
end

%% ============================================================
%  核心函数 2: 从源码提取系数
% ============================================================
function [C, inds] = extract_gsw_coeffs()
    % 定位文件
    file_path = which('gsw_specvol.m');
    if isempty(file_path)
        % 尝试旧版本或其他常见文件名
        file_path = which('gsw_specvol_t_exact.m');
    end
    
    if isempty(file_path)
        error('未找到 gsw_specvol.m 文件。请确认 GSW 已安装。');
    end
    
    % 读取文件内容
    fid = fopen(file_path, 'r');
    if fid == -1
        error('无法打开文件: %s', file_path);
    end
    file_content = fread(fid, '*char')';
    fclose(fid);
    
    % 正则表达式提取
    % 匹配模式: g000 = 1.234; 或 v000 = 1.234;
    % 同时也提取三个数字作为指数 i, j, k
    pattern = '([vg])(\d)(\d)(\d)\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?);';
    tokens = regexp(file_content, pattern, 'tokens');
    
    num_coeffs = length(tokens);
    if num_coeffs < 75
        error('未能提取到足够的系数 (仅找到 %d 个, 预期 75 个)。可能是文件格式不同。', num_coeffs);
    end
    
    % 转换为矩阵
    C = zeros(num_coeffs, 1);
    inds = zeros(num_coeffs, 3);
    
    for k = 1:num_coeffs
        tk = tokens{k};
        % tk{1} 是变量前缀 'v' 或 'g'
        % tk{2}, tk{3}, tk{4} 是指数 i, j, k
        % tk{5} 是数值字符串
        
        i = str2double(tk{2});
        j = str2double(tk{3});
        k_exp = str2double(tk{4}); % 避免与循环变量 k 冲突
        val = str2double(tk{5});
        
        C(k) = val;
        inds(k, :) = [i, j, k_exp];
    end
end