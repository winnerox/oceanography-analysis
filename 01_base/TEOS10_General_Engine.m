classdef TEOS10_General_Engine < handle
    % TEOS10_Optimized_Engine (Dynamic & Vectorized & Mixed Derivs)
    % -----------------------------------------------------------
    % [终极升级版] 支持任意阶数 (1-8阶) 的所有 纯项 + 混合项 同步计算！
    % -----------------------------------------------------------
    
    properties (Access = private)
        CoeffsMatrix % [N x 4]
        BinomialTable
        I_vec, J_vec, K_vec, Val_vec
        SCALE_S = 0.0248826675584615;
        SCALE_T = 0.025;
        SCALE_P = 1e-4;
    end
    
    methods
        %% 1. 构造函数
        function obj = TEOS10_General_Engine(max_possible_order)
            if nargin < 1, max_possible_order = 10; end 
            obj.load_gsw_coefficients();
            obj.I_vec = obj.CoeffsMatrix(:, 1);
            obj.J_vec = obj.CoeffsMatrix(:, 2);
            obj.K_vec = obj.CoeffsMatrix(:, 3);
            obj.Val_vec = obj.CoeffsMatrix(:, 4);
            obj.BinomialTable = zeros(max_possible_order+1, max_possible_order+1);
            for n = 0:max_possible_order
                for k = 0:n
                    obj.BinomialTable(n+1, k+1) = nchoosek(n, k);
                end
            end
        end
        
        %% [NEW] 全面升级：一次性计算任意阶所有的偏导数 (纯项 + 混合项)
        function [Derivs_T, Derivs_S, Derivs_Cross] = calculate_all_derivatives(obj, SA, CT, p, max_order)
            num_pts = numel(SA);
            orig_size = size(SA);
            
            SA = SA(:); CT = CT(:); p = p(:);
            
            X_S = SA + 24; X_T = CT; X_P = p;
            Iv = obj.I_vec; Jv = obj.J_vec; Kv = obj.K_vec; Val = obj.Val_vec;
            Ss = obj.SCALE_S; St = obj.SCALE_T; Sp = obj.SCALE_P;
            
            S_base = (Ss * X_S); P_base = (Sp * X_P); T_base = (St * X_T);
            P_mat = P_base .^ (Kv');
            
            % 1. 计算 v (比容) 的各阶导数网格
            v_grid = cell(max_order+1, max_order+1); 
            
            for ns = 0:max_order
                for nt = 0:max_order
                    if ns + nt > max_order, continue; end
                    
                    if ns == 0
                        mult_S = ones(size(Iv')); 
                        S_pow = Iv' ./ 2;
                    else
                        mult_S = ones(size(Iv'));
                        exp_S = Iv' ./ 2;
                        for i = 0:ns-1, mult_S = mult_S .* (exp_S - i); end
                        mult_S = mult_S .* (Ss^ns);
                        S_pow = exp_S - ns;
                    end
                    
                    if nt == 0
                        mult_T = ones(size(Jv'));
                        T_pow = Jv';
                    else
                        mult_T = ones(size(Jv'));
                        exp_T = Jv';
                        for i = 0:nt-1, mult_T = mult_T .* (exp_T - i); end
                        mult_T = mult_T .* (St^nt);
                        T_pow = exp_T - nt;
                    end
                    
                    valid_mask = (Jv' >= nt);
                    total_mult = mult_T .* mult_S .* valid_mask;
                    
                    S_term = S_base .^ S_pow;
                    T_term = T_base .^ T_pow;
                    T_term(:, ~valid_mask) = 0;
                    
                    effective_coeffs = Val .* total_mult'; 
                    v_grid{nt+1, ns+1} = (S_term .* P_mat .* T_term) * effective_coeffs;
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
        
        %% ========================================================
        %% [兼容保留原有旧方法，防止报错]
        %% ========================================================
        function out = calculate_all_orders(obj, SA, CT, p, max_order)
            if nargin < 5, max_order = 8; end
            out = struct();
            drho_dT = obj.calc_leibniz_vectorized(SA, CT, p, max_order, 'CT');
            drho_dS = obj.calc_leibniz_vectorized(SA, CT, p, max_order, 'SA');
            for n = 1:max_order
                out.(sprintf('d%d_T', n)) = drho_dT(:, n+1); 
                out.(sprintf('d%d_S', n)) = drho_dS(:, n+1);
            end
        end

        function out = calculate_rho_mixed_3rd_order(obj, SA, CT, p)
            out = struct();
            v_grid = obj.eval_v_mixed_grid_upto3(SA, CT, p);
            rho_grid = obj.calc_rho_from_v_mixed(v_grid);
            out.dS = rho_grid.d10; out.dT = rho_grid.d01;
            out.dSS = rho_grid.d20; out.dST = rho_grid.d11; out.dTT = rho_grid.d02;
            out.dSSS = rho_grid.d30; out.dSST = rho_grid.d21; out.dSTT = rho_grid.d12; out.dTTT = rho_grid.d03;
        end
    end
    
    methods (Access = private)
        function v_grid = eval_v_mixed_grid_upto3(obj, SA, CT, p)
            num_pts = numel(SA); SA = SA(:); CT = CT(:); p = p(:); v_grid = struct();
            X_S = SA + 24; X_T = CT; X_P = p; Iv = obj.I_vec; Jv = obj.J_vec; Kv = obj.K_vec; Val = obj.Val_vec;
            Ss = obj.SCALE_S; St = obj.SCALE_T; Sp = obj.SCALE_P;
            S_base = (Ss * X_S); S_mat_0 = S_base .^ (Iv' ./ 2); P_base = (Sp * X_P); P_mat = P_base .^ (Kv'); T_base = (St * X_T); T_mat_0 = T_base .^ (Jv');
            pairs = [0,0; 1,0; 0,1; 2,0; 1,1; 0,2; 3,0; 2,1; 1,2; 0,3];
            for k = 1:size(pairs,1)
                ns = pairs(k,1); nt = pairs(k,2);
                if ns == 0, mult_S = ones(size(Iv')); S_pow = Iv' ./ 2; else, mult_S = ones(size(Iv')); exp_S = Iv' ./ 2; for i = 0:ns-1, mult_S = mult_S .* (exp_S - i); end, mult_S = mult_S .* (Ss^ns); S_pow = exp_S - ns; end
                if nt == 0, mult_T = ones(size(Jv')); T_pow = Jv'; else, mult_T = ones(size(Jv')); exp_T = Jv'; for i = 0:nt-1, mult_T = mult_T .* (exp_T - i); end, mult_T = mult_T .* (St^nt); T_pow = exp_T - nt; end
                valid_mask = (Jv' >= nt); total_mult = mult_T .* mult_S .* valid_mask; S_term = S_base .^ S_pow; T_term = T_base .^ T_pow; T_term(:, ~valid_mask) = 0;
                effective_coeffs = Val .* total_mult'; term_val = (S_term .* P_mat .* T_term) * effective_coeffs; v_grid.(sprintf('d%d%d', ns, nt)) = term_val;
            end
        end
        
        function rho_grid = calc_rho_from_v_mixed(obj, v_grid)
            rho_grid = struct(); rho00 = 1 ./ v_grid.d00; rho_grid.d00 = rho00;
            pairs_order = { [1,0; 0,1], [2,0; 1,1; 0,2], [3,0; 2,1; 1,2; 0,3] };
            for k = 1:3
                current_pairs = pairs_order{k};
                for p_idx = 1:size(current_pairs,1)
                    n = current_pairs(p_idx, 1); m = current_pairs(p_idx, 2); sum_val = zeros(size(rho00));
                    for i = 0:n
                        for j = 0:m
                            if (i==n && j==m), continue; end
                            bin_n = obj.BinomialTable(n+1, i+1); bin_m = obj.BinomialTable(m+1, j+1);
                            r_ij = rho_grid.(sprintf('d%d%d', i, j)); v_rem = v_grid.(sprintf('d%d%d', n-i, m-j));
                            sum_val = sum_val + (bin_n * bin_m) * (r_ij .* v_rem);
                        end
                    end
                    rho_grid.(sprintf('d%d%d', n, m)) = -rho00 .* sum_val;
                end
            end
        end

        function rho_derivs = calc_leibniz_vectorized(obj, SA, CT, p, N, var_name)
            SA = SA(:); CT = CT(:); p = p(:); num_pts = length(SA);
            v_derivs = obj.eval_v_derivs_vectorized(SA, CT, p, N, var_name); rho_derivs = zeros(num_pts, N+1, 'like', SA);
            v0 = v_derivs(:, 1); rho_0 = 1 ./ v0; rho_derivs(:, 1) = rho_0;
            for n = 1:N
                sum_val = zeros(num_pts, 1, 'like', SA);
                for k = 0:(n-1), sum_val = sum_val + obj.BinomialTable(n+1, k+1) * (rho_derivs(:, k+1) .* v_derivs(:, (n-k)+1)); end
                rho_derivs(:, n+1) = -rho_0 .* sum_val;
            end
        end
        
        function v_vals = eval_v_derivs_vectorized(obj, SA, CT, p, N, target_var)
            num_pts = length(SA); v_vals = zeros(num_pts, N+1); X_S = SA + 24; X_T = CT; X_P = p; Iv = obj.I_vec; Jv = obj.J_vec; Kv = obj.K_vec; Val = obj.Val_vec;
            Ss = obj.SCALE_S; St = obj.SCALE_T; Sp = obj.SCALE_P; S_base = (Ss * X_S); S_mat = S_base .^ (Iv' ./ 2); P_base = (Sp * X_P); P_mat = P_base .^ (Kv'); T_base = (St * X_T); SP_mat = S_mat .* P_mat;
            if strcmp(target_var, 'CT')
                for ord = 0:N
                    mult_vec = ones(1, length(Jv)); pow_vec = Jv';
                    if ord > 0, valid_idx = (Jv >= ord); current_mult = ones(size(Jv)); for k = 0:(ord-1), current_mult = current_mult .* (Jv - k); end, mult_vec = current_mult' .* (St^ord); mult_vec(~valid_idx) = 0; pow_vec = pow_vec - ord; pow_vec(~valid_idx) = 0; end
                    T_mat = T_base .^ pow_vec; effective_coeffs = Val .* mult_vec'; v_vals(:, ord+1) = (SP_mat .* T_mat) * effective_coeffs;
                end
            elseif strcmp(target_var, 'SA')
                T_mat = T_base .^ (Jv'); TP_mat = T_mat .* P_mat;
                for ord = 0:N
                    exponent_vec = Iv' ./ 2.0; mult_vec = ones(size(exponent_vec)); pow_vec = exponent_vec;
                    if ord > 0, scale_factor = Ss^ord; current_mult = ones(size(exponent_vec)); for k = 0:(ord-1), current_mult = current_mult .* (exponent_vec - k); end, mult_vec = current_mult .* scale_factor; valid_idx = (Iv > 0); pow_vec = pow_vec - ord; end
                    S_mat_deriv = S_base .^ pow_vec; effective_coeffs = Val .* mult_vec'; v_vals(:, ord+1) = (TP_mat .* S_mat_deriv) * effective_coeffs;
                end
            end
        end
        
        function load_gsw_coefficients(obj)
            target = 'gsw_specvol.m'; fpath = which(target);
            if isempty(fpath), error('❌ 未找到 GSW 工具箱！'); end
            lines = readlines(fpath); coeffs = [];
            rgx = '^\s*v(\d)(\d)(\d)\s*=\s*([-+]?[\d\.]+(?:[eE][-+]?\d+)?);';
            for k = 1:length(lines)
                tk = regexp(lines(k), rgx, 'tokens');
                if ~isempty(tk), tokens = tk{1}; coeffs = [coeffs; str2double(tokens{1}), str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4})]; end
            end
            if isempty(coeffs), error('❌ 读取系数失败'); end
            obj.CoeffsMatrix = coeffs;
        end
    end
end