classdef gsw_rho_high_order < handle
    % TEOS10_High_Order - 任意高阶混合偏导计算器
    % =====================================================================
    % 核心升级:
    %   1. [算法重构] 废弃硬编码公式，采用 Leibniz 递归算法。
    %   2. [无限阶数] 支持 4阶, 5阶, 10阶... 只要你需要。
    %   3. [记忆化]   内部自动缓存低阶导数，避免重复计算。
    % =====================================================================
    
    properties (Access = private, Constant)
        SCALE_S = 0.0248826675584615;
        SCALE_T = 0.025;
        SCALE_P = 1e-4;
    end
    
    methods (Static)
        %% 主接口
        function result = derivative(SA, CT, p, m, n, var1, var2)
            % 计算 ∂^(m+n)ρ/∂var1^m∂var2^n
            
            if nargin < 7, error('参数不足'); end
            
            % 统一列向量
            SA = double(SA(:)); CT = double(CT(:)); p = double(p(:));
            
            % 变量解析
            [v1, v2] = gsw_rho_high_order.parse_vars(var1, var2);
            
            % 初始化记忆化存储 (Memoization)
            % 使用 Cell 数组存储中间结果: Memo{i+1, j+1} 存储 rho^(i,j)
            Memo = cell(m+1, n+1);
            
            % 开始递归求解
            result = gsw_rho_high_order.solve_recursive(m, n, v1, v2, SA, CT, p, Memo);
        end
    end
    
    methods (Static, Access = private)
        %% 核心递归求解器 (Leibniz Algorithm)
        function [res, Memo] = solve_recursive(m, n, v1, v2, SA, CT, p, Memo)
            % 检查缓存
            if ~isempty(Memo{m+1, n+1})
                res = Memo{m+1, n+1};
                return;
            end
            
            % 1. 获取 v 的当前阶导数 v^(m,n)
            v_mn = gsw_rho_high_order.calc_v_mixed(SA, CT, p, m, n, v1, v2);
            
            % 2. 递归基 (Base Case): 0阶导数
            if m == 0 && n == 0
                res = 1 ./ v_mn; % rho = 1/v
                Memo{1, 1} = res;
                return;
            end
            
            % 3. 递归步骤: Leibniz 公式
            % rho^(m,n) = -1/v * sum( C(m,i)C(n,j) * rho^(i,j) * v^(m-i,n-j) )
            
            % 先获取 v0 (即 v^(0,0)) 用于除法，为了效率可以直接调缓存
            % 这里为了代码清晰，重新获取一次(由于 persistent 缓存，开销极小)
            if m==0 && n==0
                v0 = v_mn;
            else
                v0 = gsw_rho_high_order.calc_v_mixed(SA, CT, p, 0, 0, v1, v2);
            end
            
            sum_val = zeros(size(SA));
            
            for i = 0:m
                for j = 0:n
                    % 跳过最高阶项 (自己)
                    if i == m && j == n
                        continue;
                    end
                    
                    % 计算二项式系数
                    bin_coeff = nchoosek(m, i) * nchoosek(n, j);
                    
                    % A. 递归获取 rho^(i,j)
                    [rho_ij, Memo] = gsw_rho_high_order.solve_recursive(i, j, v1, v2, SA, CT, p, Memo);
                    
                    % B. 获取 v^(m-i, n-j)
                    v_remain = gsw_rho_high_order.calc_v_mixed(SA, CT, p, m-i, n-j, v1, v2);
                    
                    % 累加
                    sum_val = sum_val + (bin_coeff * rho_ij .* v_remain);
                end
            end
            
            % 最终计算
            res = - (1 ./ v0) .* sum_val;
            
            % 存入缓存
            Memo{m+1, n+1} = res;
        end

        %% 变量解析
        function [v1, v2] = parse_vars(var1, var2)
            map = containers.Map({'S','T','P','SA','CT'}, {'S','T','P','S','T'});
            v1 = map(upper(char(var1))); % 强转char防string
            v2 = map(upper(char(var2)));
        end
        
        %% 计算 v 的混合偏导 (底层多项式引擎)
        function v_mixed = calc_v_mixed(SA, CT, p, m, n, var1, var2)
            cache = gsw_rho_high_order.get_cache();
            
            S_base = gsw_rho_high_order.SCALE_S * (SA + 24);
            T_base = gsw_rho_high_order.SCALE_T * CT;
            P_base = gsw_rho_high_order.SCALE_P * p;
            
            % 确定求导阶数
            ord_s=0; ord_t=0; ord_p=0;
            if var1=='S', ord_s=ord_s+m; elseif var1=='T', ord_t=ord_t+m; elseif var1=='P', ord_p=ord_p+m; end
            if var2=='S', ord_s=ord_s+n; elseif var2=='T', ord_t=ord_t+n; elseif var2=='P', ord_p=ord_p+n; end
            
            key = sprintf('v_%d%d%d', ord_s, ord_t, ord_p);
            
            if ~isfield(cache, key)
                cache.(key) = gsw_rho_high_order.prepare_coeffs(ord_s, ord_t, ord_p);
                gsw_rho_high_order.update_cache(cache);
            end
            
            % 向量化计算
            S_struct = cache.(key);
            idx = S_struct.valid_idx;
            if ~any(idx)
                v_mixed = zeros(size(SA));
            else
                M_S = S_base .^ S_struct.pow_s;
                M_T = T_base .^ S_struct.pow_t;
                M_P = P_base .^ S_struct.pow_p;
                v_mixed = (M_S .* M_T .* M_P) * S_struct.coeffs';
            end
        end
        
        %% 系数预计算 (与之前相同)
        function coeff_struct = prepare_coeffs(ord_s, ord_t, ord_p)
            [Iv, Jv, Kv, Gv] = gsw_rho_high_order.load_gsw_coeffs();
            Ss = gsw_rho_high_order.SCALE_S; 
            St = gsw_rho_high_order.SCALE_T; 
            Sp = gsw_rho_high_order.SCALE_P;
            
            % S 部分
            c_s = ones(size(Iv)); p_s = Iv/2;
            if ord_s>0, for k=0:ord_s-1, c_s=c_s.*(p_s-k); end; c_s=c_s*Ss^ord_s; p_s=p_s-ord_s; end
            
            % T 部分
            c_t = ones(size(Jv)); p_t = Jv;
            if ord_t>0, for k=0:ord_t-1, c_t=c_t.*(p_t-k); end; c_t=c_t*St^ord_t; p_t=p_t-ord_t; end
            
            % P 部分
            c_p = ones(size(Kv)); p_p = Kv;
            if ord_p>0, for k=0:ord_p-1, c_p=c_p.*(p_p-k); end; c_p=c_p*Sp^ord_p; p_p=p_p-ord_p; end
            
            total_c = Gv .* c_s .* c_t .* c_p;
            valid = (total_c ~= 0);
            
            coeff_struct.valid_idx = valid;
            coeff_struct.coeffs = total_c(valid)';
            coeff_struct.pow_s = p_s(valid)';
            coeff_struct.pow_t = p_t(valid)';
            coeff_struct.pow_p = p_p(valid)';
        end
        
        %% 缓存与加载 (Helper Methods)
        function cache = get_cache()
            persistent Cache; if isempty(Cache), Cache = struct(); end; cache = Cache;
        end
        function update_cache(new_cache)
            persistent Cache; Cache = new_cache;
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