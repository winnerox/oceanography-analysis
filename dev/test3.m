function teos10_final_test()
%% =========================================================================
%  TEOS-10 密度全阶偏导数计算器 (增强输出版)
%  -------------------------------------------------------------------------
%  核心功能:
%    1. [逆向工程] 自动抓取 gsw_specvol.m 系数
%    2. [智能适配] 自动识别并应用“根号缩放” (适配 GSW v3.06)
%    3. [符号求导] 计算三阶解析解
%    4. [原生输出] 显式打印 GSW 工具箱的计算结果用于核对
% =========================================================================
    clear; clc;

    % 1. 设置测试点
    SA_val = 35;       % 绝对盐度 [g/kg]
    CT_val = 20;       % 保守温度 [deg C]
    p_val  = 1000;     % 压力 [dbar]
    
    fprintf('============================================================\n');
    fprintf('  TEOS-10 密度偏导数深度验证\n');
    fprintf('  测试点: SA=%.2f, CT=%.2f, p=%.1f\n', SA_val, CT_val, p_val);
    fprintf('============================================================\n');

    %% 2. 提取系数 (Reverse Engineering)
    target_file = 'gsw_specvol.m';
    fpath = which(target_file);
    
    if isempty(fpath)
        error('错误: 未找到 %s。请确保 GSW 工具箱已安装。', target_file);
    end
    fprintf('[1/5] 读取核心文件: %s\n', fpath);
    
    lines = readlines(fpath);
    C_map = containers.Map();
    
    % 正则提取 v000 = ... 格式的系数
    for k = 1:length(lines)
        str = strtrim(lines(k));
        token = regexp(str, '^(v\d{3})\s*=\s*([-+]?[\d\.]+[eE][-+]?\d+);', 'tokens');
        if ~isempty(token)
            key = token{1}{1}; 
            val = str2double(token{1}{2});
            C_map(key) = val;
        end
    end
    
    if C_map.Count < 75
        error('提取失败: 仅找到 %d 个系数 (预期 75 个)。', C_map.Count);
    end
    fprintf('      成功提取 %d 个系数。\n', C_map.Count);

    %% 3. 确定缩放模式 (Scaling Logic)
    fprintf('[2/5] 校验缩放模式...\n');
    
    t_fac = 0.025;   % 1/40
    p_fac = 1e-4;    % 1/10000
    s_fac = 0.0248826675584615; % 1/40.188...
    
    rho_true = gsw_rho(SA_val, CT_val, p_val);
    v_true   = 1 / rho_true;
    
    % 测试根号缩放 (适配 v3.06)
    s_root = sqrt((SA_val + 24) * s_fac);
    v_root = calc_poly_val(C_map, s_root, CT_val*t_fac, p_val*p_fac);
    err_root = abs(v_root - v_true)/v_true;
    
    is_root_scaling = false;
    if err_root < 1e-9
        fprintf('      [匹配成功] 锁定【根号缩放】模式 (误差: %.2e)\n', err_root);
        is_root_scaling = true;
    else
        % 测试线性缩放
        s_lin = (SA_val + 24) * s_fac;
        v_lin = calc_poly_val(C_map, s_lin, CT_val*t_fac, p_val*p_fac);
        err_lin = abs(v_lin - v_true)/v_true;
        
        if err_lin < 1e-9
            fprintf('      [匹配成功] 锁定【线性缩放】模式 (误差: %.2e)\n', err_lin);
            is_root_scaling = false;
        else
            error('无法匹配任何缩放模式！请检查 GSW 版本。');
        end
    end

    %% 4. 符号求导 (Symbolic Differentiation)
    fprintf('[3/5] 执行符号解析求导...\n');
    syms SA CT p real
    
    if is_root_scaling
        s_sym = sqrt((SA + 24) * s_fac);
    else
        s_sym = (SA + 24) * s_fac;
    end
    t_sym = CT * t_fac;
    p_sym = p  * p_fac;
    
    % 构建比容多项式
    v_sym = 0;
    keys = C_map.keys;
    vals = C_map.values;
    for k = 1:length(keys)
        key = keys{k}; % 'v102'
        val = vals{k};
        i = str2double(key(2));
        j = str2double(key(3));
        k_idx = str2double(key(4));
        v_sym = v_sym + val * (s_sym^i) * (t_sym^j) * (p_sym^k_idx);
    end
    
    rho_sym = 1 / v_sym;
    
    % 一阶
    dr_T = diff(rho_sym, CT);
    dr_S = diff(rho_sym, SA);
    % 二阶
    dr_TT = diff(rho_sym, CT, 2);
    dr_SS = diff(rho_sym, SA, 2);
    dr_TS = diff(diff(rho_sym, CT), SA);
    % 三阶
    dr_TTT = diff(rho_sym, CT, 3);
    dr_SSS = diff(rho_sym, SA, 3);
    dr_TTS = diff(diff(rho_sym, CT, 2), SA); 
    dr_TSS = diff(diff(rho_sym, SA, 2), CT); 
    
    % 转换为函数句柄
    f_T   = matlabFunction(dr_T,  'Vars', {SA, CT, p});
    f_S   = matlabFunction(dr_S,  'Vars', {SA, CT, p});
    f_TT  = matlabFunction(dr_TT, 'Vars', {SA, CT, p});
    f_SS  = matlabFunction(dr_SS, 'Vars', {SA, CT, p});
    f_TS  = matlabFunction(dr_TS, 'Vars', {SA, CT, p});
    f_TTT = matlabFunction(dr_TTT, 'Vars', {SA, CT, p});
    f_SSS = matlabFunction(dr_SSS, 'Vars', {SA, CT, p});
    f_TTS = matlabFunction(dr_TTS, 'Vars', {SA, CT, p});
    f_TSS = matlabFunction(dr_TSS, 'Vars', {SA, CT, p});

    %% 5. 调用 GSW 工具箱原函数 (Ground Truth)
    fprintf('[4/5] 调用 GSW 原生函数获取真值...\n');
    
    % GSW 原生函数调用
    % 注意: GSW 输出顺序通常为 (SA, CT, P)
    [g_S, g_T, ~] = gsw_rho_first_derivatives(SA_val, CT_val, p_val);
    [g_SS, g_TS, g_TT, ~, ~] = gsw_rho_second_derivatives(SA_val, CT_val, p_val);

    %% 6. 最终完整输出
    fprintf('[5/5] 生成最终报告...\n');
    
    % 计算我的解析解
    my_T = f_T(SA_val, CT_val, p_val); my_S = f_S(SA_val, CT_val, p_val);
    my_TT = f_TT(SA_val, CT_val, p_val); my_SS = f_SS(SA_val, CT_val, p_val); my_TS = f_TS(SA_val, CT_val, p_val);
    my_TTT = f_TTT(SA_val, CT_val, p_val); my_SSS = f_SSS(SA_val, CT_val, p_val);
    my_TTS = f_TTS(SA_val, CT_val, p_val); my_TSS = f_TSS(SA_val, CT_val, p_val);

    fprintf('\n==========================================================================\n');
    fprintf('  TEOS-10 密度偏导数全集 (My Calc vs GSW Toolbox)\n');
    fprintf('  (单位: SI 标准单位, 相对误差需 < 1e-10)\n');
    fprintf('==========================================================================\n');
    
    % 打印一阶导数
    fprintf(' [一阶导数]\n');
    print_compare('rho_T (对温度)', my_T, g_T);
    print_compare('rho_S (对盐度)', my_S, g_S);
    
    % 打印二阶导数
    fprintf('\n [二阶导数]\n');
    print_compare('rho_TT (温度二阶)', my_TT, g_TT);
    print_compare('rho_SS (盐度二阶)', my_SS, g_SS);
    print_compare('rho_TS (混合二阶)', my_TS, g_TS);
    
    % 打印三阶导数 (仅 My Calc，因为 GSW 没有)
    fprintf('\n [三阶导数] (GSW 无此功能，以下为基于底层一致性的解析真值)\n');
    fprintf('  rho_TTT : % .12e\n', my_TTT);
    fprintf('  rho_SSS : % .12e\n', my_SSS);
    fprintf('  rho_TTS : % .12e\n', my_TTS);
    fprintf('  rho_TSS : % .12e\n', my_TSS);
    fprintf('==========================================================================\n');
end

% 辅助函数: 打印对比行
function print_compare(label, my_val, gsw_val)
    err = abs((my_val - gsw_val)/gsw_val);
    % 格式化输出: 标签 | 我的计算 | GSW真值 | 误差
    fprintf('  %-18s | My: % .8e | GSW: % .8e | Err: %.2e\n', ...
        label, my_val, gsw_val, err);
end

% 辅助函数: 多项式计算
function v = calc_poly_val(C_map, s, y, z)
    v = 0;
    keys = C_map.keys;
    vals = C_map.values;
    for k = 1:length(keys)
        key = keys{k};
        val = vals{k};
        i = str2double(key(2));
        j = str2double(key(3));
        k_idx = str2double(key(4));
        v = v + val * (s^i) * (y^j) * (z^k_idx);
    end
end