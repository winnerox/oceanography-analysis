classdef TEOS10_Density_Engine < handle
    % TEOS10_Density_Engine
    % =========================================================================
    % 纯正的【密度 (Density, rho)】高阶求导引擎
    %
    % 核心算法：
    % 1. 解析 TEOS-10 官方的比容多项式 (gsw_specvol.m) 获得 v 的高阶导数。
    % 2. [突破] 利用莱布尼茨求导法则 (v * rho = 1)，通过严格的递归方程：
    %    rho^(n) = -(1/v) * sum_{k=1}^n [ C(n,k) * v^(k) * rho^(n-k) ]
    %    精确、无损地将比容导数转换为真实的【密度高阶导数】。
    % =========================================================================
    
    properties (Access = private)
        I_vec
        J_vec
        K_vec
        Val_vec
        NumTerms
    end
    
    methods
        %% 1. 构造函数
        function obj = TEOS10_Density_Engine(~)
            obj.load_gsw_coefficients();
        end
        
        %% 2. 解析 GSW 工具箱系数
        function load_gsw_coefficients(obj)
            target = 'gsw_specvol.m';
            fpath = which(target);
            if isempty(fpath)
                error('❌ 未找到 GSW 工具箱！请确保已将 GSW 添加到 MATLAB 路径中。');
            end
            
            lines = readlines(fpath);
            rgx = '^\s*v(\d)(\d)(\d)\s*=\s*([-+]?[\d\.]+(?:[eE][-+]?\d+)?);';
            
            i_list = []; j_list = []; k_list = []; val_list = [];
            
            for idx = 1:length(lines)
                tk = regexp(lines(idx), rgx, 'tokens');
                if ~isempty(tk)
                    tokens = tk{1};
                    i_list(end+1) = str2double(tokens{1});     % SA 幂次
                    j_list(end+1) = str2double(tokens{2});     % CT 幂次
                    k_list(end+1) = str2double(tokens{3});     % P  幂次
                    val_list(end+1) = str2double(tokens{4});   % 系数
                end
            end
            
            if isempty(val_list)
                error('❌ 未能从 gsw_specvol.m 提取任何系数。');
            end
            
            obj.I_vec = i_list(:);
            obj.J_vec = j_list(:);
            obj.K_vec = k_list(:);
            obj.Val_vec = val_list(:);
            obj.NumTerms = length(val_list);
            
            fprintf('   [引擎] 成功加载多项式系数。已开启【密度(rho)导数转换模式】。\n');
        end
        
        %% 3. 计算【密度(rho)】对温度的全阶偏导数
        function rho_derivs = calculate_T_derivs(obj, SA, CT, p, max_order)
            N = length(SA);
            v_derivs = zeros(N, max_order + 1, 'like', SA);
            
            % Step 1: 先算出比容 (v) 的各阶偏导数
            for ord = 0:max_order
                v_derivs(:, ord+1) = obj.evaluate_deriv(SA, CT, p, ord, 0);
            end
            
            % Step 2: 莱布尼茨法则精确转换为密度 (rho) 的各阶偏导数
            rho_derivs = zeros(N, max_order + 1, 'like', SA);
            rho_derivs(:, 1) = 1 ./ v_derivs(:, 1); % 0阶：rho = 1/v
            
            for n = 1:max_order
                sum_term = zeros(N, 1, 'like', SA);
                for k = 1:n
                    nCk = nchoosek(n, k);
                    sum_term = sum_term + nCk .* v_derivs(:, k+1) .* rho_derivs(:, n-k+1);
                end
                rho_derivs(:, n+1) = - sum_term ./ v_derivs(:, 1);
            end
        end
        
        %% 4. 计算【密度(rho)】对盐度的全阶偏导数
        function rho_derivs = calculate_S_derivs(obj, SA, CT, p, max_order)
            N = length(SA);
            v_derivs = zeros(N, max_order + 1, 'like', SA);
            
            % Step 1: 先算出比容 (v) 的各阶偏导数
            for ord = 0:max_order
                v_derivs(:, ord+1) = obj.evaluate_deriv(SA, CT, p, 0, ord);
            end
            
            % Step 2: 莱布尼茨法则精确转换为密度 (rho) 的各阶偏导数
            rho_derivs = zeros(N, max_order + 1, 'like', SA);
            rho_derivs(:, 1) = 1 ./ v_derivs(:, 1); % 0阶：rho = 1/v
            
            for n = 1:max_order
                sum_term = zeros(N, 1, 'like', SA);
                for k = 1:n
                    nCk = nchoosek(n, k);
                    sum_term = sum_term + nCk .* v_derivs(:, k+1) .* rho_derivs(:, n-k+1);
                end
                rho_derivs(:, n+1) = - sum_term ./ v_derivs(:, 1);
            end
        end
        
    end
    
    methods (Access = private)
        %% 5. 核心多项式向量化求导极速算子 (针对比容 v)
        function result = evaluate_deriv(obj, SA, CT, p, nT, nS)
            N = length(SA);
            result = zeros(N, 1, 'like', SA);
            
            for c = 1:obj.NumTerms
                i = obj.I_vec(c); j = obj.J_vec(c); k = obj.K_vec(c); val = obj.Val_vec(c);
                if j >= nT && i >= nS
                    coeff = val;
                    for step = 0:(nT-1), coeff = coeff * (j - step); end
                    for step = 0:(nS-1), coeff = coeff * (i - step); end
                    term = coeff .* (SA .^ (i - nS)) .* (CT .^ (j - nT)) .* (p .^ k);
                    result = result + term;
                end
            end
        end
    end
end