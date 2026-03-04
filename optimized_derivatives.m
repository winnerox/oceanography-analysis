%% 优化版导数计算 - 减少重复计算
function [Derivs_T, Derivs_S, Derivs_Cross] = calculate_all_derivatives_optimized(obj, SA, CT, p, max_order)
    % 优化版本：预计算所有幂次，避免重复计算
    
    num_pts = numel(SA);
    orig_size = size(SA);
    
    SA = SA(:); CT = CT(:); p = p(:);
    
    X_S = SA + 24; X_T = CT; X_P = p;
    Iv = obj.I_vec; Jv = obj.J_vec; Kv = obj.K_vec; Val = obj.Val_vec;
    Ss = obj.SCALE_S; St = obj.SCALE_T; Sp = obj.SCALE_P;
    
    % 预计算所有需要的幂次
    S_base = (Ss * X_S);
    P_base = (Sp * X_P);
    T_base = (St * X_T);
    
    % 预计算 S 的所有幂次 (0 到 max_order)
    S_powers = cell(max_order+1, 1);
    for ns = 0:max_order
        S_pow_vec = Iv' ./ 2 - ns;
        S_powers{ns+1} = S_base .^ S_pow_vec;
    end
    
    % 预计算 T 的所有幂次 (0 到 max_order)
    T_powers = cell(max_order+1, 1);
    for nt = 0:max_order
        T_pow_vec = Jv' - nt;
        valid_mask = (Jv' >= nt);
        T_term = T_base .^ T_pow_vec;
        T_term(:, ~valid_mask) = 0;
        T_powers{nt+1} = T_term;
    end
    
    % 预计算 P 的幂次
    P_mat = P_base .^ (Kv');
    
    % 预计算系数乘子
    S_multipliers = cell(max_order+1, 1);
    for ns = 0:max_order
        if ns == 0
            S_multipliers{ns+1} = ones(size(Iv'));
        else
            mult_S = ones(size(Iv'));
            exp_S = Iv' ./ 2;
            for i = 0:ns-1
                mult_S = mult_S .* (exp_S - i);
            end
            S_multipliers{ns+1} = mult_S .* (Ss^ns);
        end
    end
    
    T_multipliers = cell(max_order+1, 1);
    for nt = 0:max_order
        if nt == 0
            T_multipliers{nt+1} = ones(size(Jv'));
        else
            mult_T = ones(size(Jv'));
            exp_T = Jv';
            for i = 0:nt-1
                mult_T = mult_T .* (exp_T - i);
            end
            T_multipliers{nt+1} = mult_T .* (St^nt);
        end
    end
    
    % 1. 计算 v (比容) 的各阶导数网格
    v_grid = cell(max_order+1, max_order+1);
    
    for ns = 0:max_order
        for nt = 0:max_order
            if ns + nt > max_order, continue; end
            
            % 使用预计算的结果
            S_term = S_powers{ns+1};
            T_term = T_powers{nt+1};
            mult_S = S_multipliers{ns+1};
            mult_T = T_multipliers{nt+1};
            
            valid_mask = (Jv' >= nt);
            total_mult = mult_T .* mult_S .* valid_mask;
            
            effective_coeffs = Val .* total_mult';
            term_val = (S_term .* P_mat .* T_term) * effective_coeffs;
            
            v_grid{nt+1, ns+1} = term_val;
        end
    end
    
    % 2. 广义莱布尼茨反演求 rho 密度导数
    rho_grid = cell(max_order+1, max_order+1);
    rho00 = 1 ./ v_grid{1, 1};
    rho_grid{1, 1} = rho00;
    
    for total_order = 1:max_order
        for nt = 0:total_order
            ns = total_order - nt;
            
            sum_val = zeros(num_pts, 1, 'like', SA);
            for i = 0:nt
                for j = 0:ns
                    if (i == nt && j == ns), continue; end
                    bin_nt = obj.BinomialTable(nt+1, i+1);
                    bin_ns = obj.BinomialTable(ns+1, j+1);
                    
                    r_ij = rho_grid{i+1, j+1};
                    v_rem = v_grid{nt-i+1, ns-j+1};
                    sum_val = sum_val + (bin_nt * bin_ns) * (r_ij .* v_rem);
                end
            end
            rho_grid{nt+1, ns+1} = -rho00 .* sum_val;
        end
    end
    
    % 3. 提取所需的导数并恢复原始维度
    Derivs_T = cell(max_order, 1);
    Derivs_S = cell(max_order, 1);
    Derivs_Cross = cell(max_order, max_order);
    
    for n = 1:max_order
        Derivs_T{n} = reshape(rho_grid{n+1, 1}, orig_size);
        Derivs_S{n} = reshape(rho_grid{1, n+1}, orig_size);
        for nt = 1:(n-1)
            ns = n - nt;
            Derivs_Cross{nt, ns} = reshape(rho_grid{nt+1, ns+1}, orig_size);
        end
    end
end