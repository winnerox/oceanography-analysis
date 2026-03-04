classdef TEOS10_Lower_Derivatives < handle
    % TEOS10_Lower_Derivatives 
    % 用于计算 TEOS-10 密度的一阶和二阶偏导数，并支持与 GSW 工具箱进行精度对齐
    
    properties (SetAccess = private)
        IsReady = false;           % 初始化状态
        ScalingMode = '';          % 自动检测到的缩放模式
        CoeffsMap = [];            % 系数缓存
        FuncHandles = struct();    % 编译后的函数句柄 (计算速度极快)
    end
    
    methods
        %% 构造函数
        function obj = TEOS10_Lower_Derivatives()
            try
                % 1. 提取系数 (与三阶导数逻辑一致)
                obj.CoeffsMap = obj.extract_coefficients();
                
                % 2. 构建符号模型 (核心：处理根号缩放)
                [rho_sym, vars] = obj.build_symbolic_model();
                
                % 3. 编译一阶和二阶导数公式
                obj.compile_derivatives(rho_sym, vars);
                
                obj.IsReady = true;
                fprintf('TEOS-10 低阶导数验证引擎已就绪 (%s)。\n', obj.ScalingMode);
            catch ME
                error('引擎初始化失败: %s', ME.message);
            end
        end
        
        %% 计算方法: 获取一阶和二阶导数
        function out = calculate(obj, SA, CT, p)
            if ~obj.IsReady, error('引擎未初始化'); end
            
            % 维度适配
            if isscalar(p) && ~isscalar(SA)
                p = repmat(p, size(SA));
            end
            
            % 1. 一阶导数
            out.rho_SA = obj.FuncHandles.rho_SA(SA, CT, p);
            out.rho_CT = obj.FuncHandles.rho_CT(SA, CT, p);
            
            % 2. 二阶导数
            out.rho_SA_SA = obj.FuncHandles.rho_SA_SA(SA, CT, p);
            out.rho_CT_CT = obj.FuncHandles.rho_CT_CT(SA, CT, p);
            out.rho_SA_CT = obj.FuncHandles.rho_SA_CT(SA, CT, p);
        end
        
        %% 核心功能: 与 GSW 工具箱进行对比验证
        function compare_with_gsw(obj, SA, CT, p)
            % 运行此方法前，请确保安装了 GSW 工具箱
            
            fprintf('\n========== TEOS-10 导数精度对比报告 ==========\n');
            fprintf('测试条件: SA=%.2f g/kg, CT=%.2f degC, p=%.0f dbar\n', ...
                mean(SA(:)), mean(CT(:)), mean(p(:)));
            fprintf('------------------------------------------------------\n');
            fprintf('%-10s | %-15s | %-15s | %-12s\n', '导数项', '本引擎结果', 'GSW工具箱结果', '相对误差');
            fprintf('------------------------------------------------------\n');
            
            % 1. 获取本引擎结果
            my_res = obj.calculate(SA, CT, p);
            
            % 2. 获取 GSW 结果
            % GSW 一阶: [rho_SA, rho_CT, rho_P]
            [gsw_SA, gsw_CT, ~] = gsw_rho_first_derivatives(SA, CT, p);
            
            % GSW 二阶: [rho_SA_SA, rho_SA_CT, rho_CT_CT, rho_SA_P, rho_CT_P]
            [gsw_SA_SA, gsw_SA_CT, gsw_CT_CT, ~, ~] = gsw_rho_second_derivatives(SA, CT, p);
            
            % 3. 逐项对比
            obj.print_diff('rho_SA',    my_res.rho_SA,    gsw_SA);
            obj.print_diff('rho_CT',    my_res.rho_CT,    gsw_CT);
            obj.print_diff('rho_SA_SA', my_res.rho_SA_SA, gsw_SA_SA);
            obj.print_diff('rho_CT_CT', my_res.rho_CT_CT, gsw_CT_CT);
            obj.print_diff('rho_SA_CT', my_res.rho_SA_CT, gsw_SA_CT);
            
            fprintf('------------------------------------------------------\n');
        end
    end
    
    methods (Access = private)
        %% 辅助: 打印差异
        function print_diff(~, name, val_mine, val_gsw)
            % 选取第一个点做展示
            v1 = val_mine(1);
            v2 = val_gsw(1);
            % 计算相对误差 (加 eps 防止除零)
            rel_err = abs(v1 - v2) / (abs(v2) + eps);
            
            % 颜色标记: 绿色优秀，红色警告
            if rel_err < 1e-12
                status = '[完美]';
            elseif rel_err < 1e-9
                status = '[合格]';
            else
                status = '[警告!]';
            end
            
            fprintf('%-10s | %15.8e | %15.8e | %10.2e %s\n', ...
                name, v1, v2, rel_err, status);
        end
        
        %% 内部: 提取系数 (复用)
        function map = extract_coefficients(~)
            target = 'gsw_specvol.m';
            fpath = which(target);
            if isempty(fpath), error('未找到 GSW 工具箱 (%s)', target); end
            lines = readlines(fpath);
            map = containers.Map();
            for k = 1:length(lines)
                str = strtrim(lines(k));
                tk = regexp(str, '^(v\d{3})\s*=\s*([-+]?[\d\.]+[eE][-+]?\d+);', 'tokens');
                if ~isempty(tk), map(tk{1}{1}) = str2double(tk{1}{2}); end
            end
        end
        
        %% 内部: 构建模型 (复用，确保 Scaling 一致)
        function [rho_sym, vars] = build_symbolic_model(obj)
            syms SA CT p real
            s_fac = 0.0248826675584615; t_fac = 0.025; p_fac = 1e-4;
            
            % 自动检测缩放模式
            sa_t = 35; ct_t = 20; p_t = 1000;
            v_gsw = 1/gsw_rho(sa_t, ct_t, p_t);
            calc_v = @(is_root) obj.calc_numeric(is_root, sa_t, ct_t, p_t);
            
            if abs(calc_v(true) - v_gsw)/v_gsw < 1e-10
                obj.ScalingMode = 'Root_Scaling (Sqrt)';
                s_term = sqrt((SA + 24) * s_fac);
            else
                obj.ScalingMode = 'Linear_Scaling';
                s_term = (SA + 24) * s_fac;
            end
            
            % 构建多项式 v
            v_sym = 0;
            keys = obj.CoeffsMap.keys; vals = obj.CoeffsMap.values;
            t_term = CT * t_fac; p_term = p * p_fac;
            for k = 1:length(keys)
                key = keys{k}; val = vals{k};
                i = str2double(key(2)); j = str2double(key(3)); k_idx = str2double(key(4));
                v_sym = v_sym + val * (s_term^i) * (t_term^j) * (p_term^k_idx);
            end
            
            rho_sym = 1 / v_sym; % 密度 = 1/比容
            vars = {SA, CT, p};
        end
        
        %% 内部: 编译导数 (针对一阶和二阶)
        function compile_derivatives(obj, rho_sym, vars)
            SA = vars{1}; CT = vars{2};
            
            % --- 1. 一阶导数 ---
            % d(rho)/dSA
            obj.FuncHandles.rho_SA = matlabFunction(diff(rho_sym, SA), 'Vars', vars);
            % d(rho)/dCT
            obj.FuncHandles.rho_CT = matlabFunction(diff(rho_sym, CT), 'Vars', vars);
            
            % --- 2. 二阶导数 ---
            % d^2(rho)/dSA^2
            obj.FuncHandles.rho_SA_SA = matlabFunction(diff(rho_sym, SA, 2), 'Vars', vars);
            % d^2(rho)/dCT^2
            obj.FuncHandles.rho_CT_CT = matlabFunction(diff(rho_sym, CT, 2), 'Vars', vars);
            % d^2(rho)/dSA dCT (混合偏导)
            obj.FuncHandles.rho_SA_CT = matlabFunction(diff(diff(rho_sym, SA), CT), 'Vars', vars);
        end
        
        %% 辅助: 数值计算用于检测缩放
        function v = calc_numeric(obj, is_root, SA, CT, p)
            s_fac = 0.0248826675584615; t_fac = 0.025; p_fac = 1e-4;
            if is_root, s = sqrt((SA + 24) * s_fac); else, s = (SA + 24) * s_fac; end
            y = CT * t_fac; z = p * p_fac;
            v = 0;
            keys = obj.CoeffsMap.keys; vals = obj.CoeffsMap.values;
            for k = 1:length(keys)
                key = keys{k}; val = vals{k};
                i = str2double(key(2)); j = str2double(key(3)); k_idx = str2double(key(4));
                v = v + val * (s^i) * (y^j) * (z^k_idx);
            end
        end
    end
end