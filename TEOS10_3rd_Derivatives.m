classdef TEOS10_3rd_Derivatives < handle
    % TEOS10_3rd_Derivatives  TEOS-10 密度三阶偏导数专用计算引擎
    %
    % 核心功能:
    %   1. 自动适配 GSW v3.06+ (自动识别根号缩放)
    %   2. 仅输出密度对温度(CT)和盐度(SA)的各项三阶偏导数
    %
    % 输出结构体包含:
    %   .rho_TTT  (d^3 rho / dCT^3)
    %   .rho_SSS  (d^3 rho / dSA^3)
    %   .rho_TTS  (d^3 rho / dCT^2 dSA)
    %   .rho_TSS  (d^3 rho / dCT dSA^2)
    
    properties (SetAccess = private)
        IsReady = false;           % 初始化状态
        ScalingMode = '';          % 缩放模式记录
        CoeffsMap = [];            % 系数容器
        FuncHandles = struct();    % 函数句柄缓存
    end

    
    
    methods
        %% 构造函数 (初始化)
        function obj = TEOS10_3rd_Derivatives()
            try
                % 1. 提取系数
                obj.CoeffsMap = obj.extract_coefficients();
                
                % 2. 确定缩放并构建符号模型
                [rho_sym, vars] = obj.build_symbolic_model();
                
                % 3. 编译导数公式 (含自检用的低阶项)
                obj.compile_derivatives(rho_sym, vars);
                
                % 4. 启动自检 (确保环境正确)
                obj.internal_verify();
                
                obj.IsReady = true;
                fprintf('TEOS-10 三阶导数引擎已就绪 (%s)。\n', obj.ScalingMode);
                
            catch ME
                error('引擎初始化失败: %s', ME.message);
            end
        end
        
        %% 核心方法: 计算三阶导数
        function out = calculate(obj, SA, CT, p)
            % calculate仅返回三阶导数
            % 输入: SA(g/kg), CT(deg C), p(dbar)
            
            if ~obj.IsReady, error('引擎未初始化'); end
            
            % 维度检查
            if ~isequal(size(SA), size(CT))
                 if isscalar(p), p = repmat(p, size(SA)); end
            end

            % 仅计算并返回三阶偏导数
            out.rho_TTT = obj.FuncHandles.rho_TTT(SA, CT, p);
            out.rho_SSS = obj.FuncHandles.rho_SSS(SA, CT, p);
            out.rho_TTS = obj.FuncHandles.rho_TTS(SA, CT, p);
            out.rho_TSS = obj.FuncHandles.rho_TSS(SA, CT, p);
        end

       
       
    end
    
    methods (Access = private)
        %% 内部: 提取系数
        function map = extract_coefficients(~)
            target = 'gsw_specvol.m';
            fpath = which(target);
            if isempty(fpath), error('未找到 %s', target); end
            
            lines = readlines(fpath);
            map = containers.Map();
            for k = 1:length(lines)
                str = strtrim(lines(k));
                tk = regexp(str, '^(v\d{3})\s*=\s*([-+]?[\d\.]+[eE][-+]?\d+);', 'tokens');
                if ~isempty(tk), map(tk{1}{1}) = str2double(tk{1}{2}); end
            end
            if map.Count < 75, error('系数提取不完整 (%d/75)', map.Count); end
        end
        
        %% 内部: 构建模型
        function [rho_sym, vars] = build_symbolic_model(obj)
            syms SA CT p real
            s_fac = 0.0248826675584615; t_fac = 0.025; p_fac = 1e-4;
            
            % 逆向检测
            sa_t = 35; ct_t = 20; p_t = 1000;
            v_true = 1/gsw_rho(sa_t, ct_t, p_t);
            
            calc_v = @(is_root) obj.calc_numeric(is_root, sa_t, ct_t, p_t);
            
            if abs(calc_v(true) - v_true)/v_true < 1e-10
                obj.ScalingMode = 'Root_Scaling';
                s_term = sqrt((SA + 24) * s_fac);
            else
                obj.ScalingMode = 'Linear_Scaling';
                s_term = (SA + 24) * s_fac;
            end
            
            % 构建多项式
            v_sym = 0;
            keys = obj.CoeffsMap.keys; vals = obj.CoeffsMap.values;
            t_term = CT * t_fac; p_term = p * p_fac;
            
            for k = 1:length(keys)
                key = keys{k}; val = vals{k};
                i = str2double(key(2)); j = str2double(key(3)); k_idx = str2double(key(4));
                v_sym = v_sym + val * (s_term^i) * (t_term^j) * (p_term^k_idx);
            end
            rho_sym = 1 / v_sym;
            vars = {SA, CT, p};
        end
        
        %% 内部: 编译导数
        function compile_derivatives(obj, rho_sym, vars)
            % 三阶导数 (输出用)
            obj.FuncHandles.rho_TTT = matlabFunction(diff(rho_sym, vars{2}, 3), 'Vars', vars);
            obj.FuncHandles.rho_SSS = matlabFunction(diff(rho_sym, vars{1}, 3), 'Vars', vars);
            obj.FuncHandles.rho_TTS = matlabFunction(diff(diff(rho_sym, vars{2}, 2), vars{1}), 'Vars', vars);
            obj.FuncHandles.rho_TSS = matlabFunction(diff(diff(rho_sym, vars{1}, 2), vars{2}), 'Vars', vars);
            
            % 一阶导数 (仅用于内部自检)
            obj.FuncHandles.rho_T = matlabFunction(diff(rho_sym, vars{2}), 'Vars', vars);
        end
        
        %% 内部: 自检
        function internal_verify(obj)
            % 隐式验证一阶导数，确保模型正确
            SA = 35; CT = 20; p = 1000;
            [~, g_T, ~] = gsw_rho_first_derivatives(SA, CT, p);
            my_T = obj.FuncHandles.rho_T(SA, CT, p);
            if abs((my_T - g_T)/g_T) > 1e-9
                warning('模型自检未通过，结果可能不准确！');
            end
        end
        
        %% 辅助: 数值计算
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