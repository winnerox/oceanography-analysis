classdef gsw_rho_mixed_derivatives < handle
    % GSW_RHO_MIXED_DERIVATIVES (Universal Compatible Version)
    % =====================================================================
    % 核心算法：迭代动态规划 (Iterative DP) - 保持不变
    % 兼容性修复：将 .^ 广播改为 bsxfun，确保 R2016b 以前版本也能运行。
    % 计算逻辑：严格保持原版数学公式，无任何改动。
    % =====================================================================
    
    properties (Access = private, Constant)
        SCALE_S = 0.0248826675584615;
        SCALE_T = 0.025;
        SCALE_P = 1e-4;
    end
    
    methods (Static)
        function result = derivative(SA, CT, p, m, n, var1, var2)
            %% 1. 预处理与内存分配
            SA = double(SA(:)); CT = double(CT(:)); p = double(p(:));
            n_points = length(SA);
            
            % 解析变量 ID
            [dim1, dim2] = gsw_rho_mixed_derivatives.parse_vars_to_id(var1, var2);
            
            % 预分配 3D 结果矩阵
            R = zeros(n_points, m+1, n+1);
            V = zeros(n_points, m+1, n+1);
            
            % 计算二项式系数
            B_m = gsw_rho_mixed_derivatives.get_binomial_coeffs(m);
            B_n = gsw_rho_mixed_derivatives.get_binomial_coeffs(n);
            
            %% 2. 批量计算所有需要的 v 的偏导数
            coeffs_map = cell(m+1, n+1);
            for i = 0:m
                for j = 0:n
                    coeffs_map{i+1, j+1} = gsw_rho_mixed_derivatives.get_poly_coeffs(i, j, dim1, dim2);
                end
            end
            
            % 基础变量缩放
            S_base = gsw_rho_mixed_derivatives.SCALE_S * (SA + 24);
            T_base = gsw_rho_mixed_derivatives.SCALE_T * CT;
            P_base = gsw_rho_mixed_derivatives.SCALE_P * p;
            
            for i = 0:m
                for j = 0:n
                    c_struct = coeffs_map{i+1, j+1};
                    if ~isempty(c_struct.coeffs)
                        % ===【兼容性修复：使用 bsxfun】===
                        % 逻辑：计算 (S^pow_s * T^pow_t * P^pow_p)
                        % 这里的数学逻辑与 .^ 完全一致，只是换了函数写法以兼容旧版 MATLAB
                        
                        % 1. 确保指数是行向量 (1 x K)
                        ps = c_struct.pow_s(:)'; 
                        pt = c_struct.pow_t(:)';
                        pp = c_struct.pow_p(:)';
                        
                        % 2. 使用底层函数 bsxfun 进行幂运算
                        % term_s(N x K) = S_base(N x 1) .^ ps(1 x K)
                        term_s = bsxfun(@power, S_base, ps); 
                        term_t = bsxfun(@power, T_base, pt);
                        term_p = bsxfun(@power, P_base, pp);
                        
                        % 3. 逐元素相乘 (N x K)
                        term_NK = term_s .* term_t .* term_p;
                        
                        % 4. 加权求和 (N x K) * (K x 1) -> (N x 1)
                        V(:, i+1, j+1) = term_NK * c_struct.coeffs;
                    end
                end
            end
            
            inv_v0 = 1 ./ V(:, 1, 1);
            
            %% 3. 动态规划求解 Rho (Core DP Loop)
            for i = 0:m
                for j = 0:n
                    if i==0 && j==0
                        R(:, 1, 1) = inv_v0;
                        continue;
                    end
                    
                    sum_val = zeros(n_points, 1);
                    
                    for i_sub = 0:i
                        for j_sub = 0:j
                            if i_sub == i && j_sub == j
                                continue;
                            end
                            coef = B_m(i+1, i_sub+1) * B_n(j+1, j_sub+1);
                            rho_prev = R(:, i_sub+1, j_sub+1);
                            v_rem = V(:, (i-i_sub)+1, (j-j_sub)+1);
                            sum_val = sum_val + coef * (rho_prev .* v_rem);
                        end
                    end
                    R(:, i+1, j+1) = - inv_v0 .* sum_val;
                end
            end
            
            result = R(:, m+1, n+1);
        end
    end
    
    methods (Static, Access = private)
        function s = get_poly_coeffs(m, n, dim1, dim2)
            ord = [0, 0, 0];
            ord(dim1) = ord(dim1) + m;
            ord(dim2) = ord(dim2) + n;
            
            cache = gsw_rho_mixed_derivatives.get_global_cache();
            num_key = ord(1)*100 + ord(2)*10 + ord(3); 
            
            if isKey(cache, num_key)
                s = cache(num_key);
            else
                s = gsw_rho_mixed_derivatives.prepare_raw_coeffs(ord(1), ord(2), ord(3));
                cache(num_key) = s;
            end
        end
        
        function B = get_binomial_coeffs(max_n)
            B = zeros(max_n+1, max_n+1);
            for i = 0:max_n
                for j = 0:i
                    B(i+1, j+1) = nchoosek(i, j);
                end
            end
        end

        function [d1, d2] = parse_vars_to_id(v1, v2)
            map = struct('S',1, 'SA',1, 'T',2, 'CT',2, 'P',3);
            d1 = map.(upper(char(v1)));
            d2 = map.(upper(char(v2)));
        end
        
        function coeff_struct = prepare_raw_coeffs(ord_s, ord_t, ord_p)
            [Iv, Jv, Kv, Gv] = gsw_rho_mixed_derivatives.load_gsw_coeffs();
            Ss = gsw_rho_mixed_derivatives.SCALE_S; 
            St = gsw_rho_mixed_derivatives.SCALE_T; 
            Sp = gsw_rho_mixed_derivatives.SCALE_P;
            
            c_s = ones(size(Iv)); p_s = Iv/2;
            if ord_s>0, for k=0:ord_s-1, c_s=c_s.*(p_s-k); end; c_s=c_s*Ss^ord_s; p_s=p_s-ord_s; end
            
            c_t = ones(size(Jv)); p_t = Jv;
            if ord_t>0, for k=0:ord_t-1, c_t=c_t.*(p_t-k); end; c_t=c_t*St^ord_t; p_t=p_t-ord_t; end
            
            c_p = ones(size(Kv)); p_p = Kv;
            if ord_p>0, for k=0:ord_p-1, c_p=c_p.*(p_p-k); end; c_p=c_p*Sp^ord_p; p_p=p_p-ord_p; end
            
            total_c = Gv .* c_s .* c_t .* c_p;
            valid = (total_c ~= 0);
            
            coeff_struct.coeffs = total_c(valid); 
            coeff_struct.pow_s = p_s(valid); 
            coeff_struct.pow_t = p_t(valid); 
            coeff_struct.pow_p = p_p(valid);
        end
        
        function c = get_global_cache()
            persistent CacheMap; if isempty(CacheMap), CacheMap = containers.Map('KeyType', 'double', 'ValueType', 'any'); end; c = CacheMap;
        end
        
        function [Iv, Jv, Kv, Gv] = load_gsw_coeffs()
            persistent I_s J_s K_s G_s
            if isempty(I_s)
                fpath = which('gsw_specvol');
                if isempty(fpath), error('需安装GSW工具箱'); end
                fid=fopen(fpath,'r'); txt=fscanf(fid,'%c'); fclose(fid);
                tk = regexp(txt, 'v(\d)(\d)(\d)\s*=\s*([-+]?[\d\.]+(?:[eE][-+]?\d+)?)', 'tokens');
                G=zeros(length(tk),4); for i=1:length(tk), G(i,:)=str2double(tk{i}); end
                I_s=G(:,1); J_s=G(:,2); K_s=G(:,3); G_s=G(:,4);
            end
            Iv=I_s; Jv=J_s; Kv=K_s; Gv=G_s;
        end
    end
end