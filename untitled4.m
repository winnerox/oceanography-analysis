function Print_GSW_Coeffs()
%% =========================================================================
%  GSW 75项多项式系数查看器
%  功能: 从 gsw_specvol.m 提取系数，按幂次 (i,j,k) 排序并打印
% =========================================================================
    clear; clc;

    % 1. 定位文件
    target_file = 'gsw_specvol.m';
    fpath = which(target_file);
    
    if isempty(fpath)
        error('未找到 gsw_specvol.m，请检查 GSW 工具箱路径。');
    end
    fprintf('正在从文件提取系数: %s\n', fpath);
    fprintf('------------------------------------------------------------\n');

    % 2. 读取并正则提取
    file_content = fileread(fpath);
    
    % 正则表达式匹配: v102 = 1.234e-5;
    % 捕获组: (i) (j) (k) (数值)
    pattern = 'v(\d)(\d)(\d)\s*=\s*([-+]?[\d\.]+[eE][-+]?\d+);';
    [tokens, ~] = regexp(file_content, pattern, 'tokens', 'match');
    
    num_coeffs = length(tokens);
    if num_coeffs == 0
        error('未提取到任何系数。文件格式可能不符。');
    end
    
    % 3. 整理数据到表格 (方便排序)
    % 预分配
    names = cell(num_coeffs, 1);
    idx_i = zeros(num_coeffs, 1);
    idx_j = zeros(num_coeffs, 1);
    idx_k = zeros(num_coeffs, 1);
    values = zeros(num_coeffs, 1);
    
    for n = 1:num_coeffs
        tk = tokens{n};
        i = str2double(tk{1});
        j = str2double(tk{2});
        k_val = str2double(tk{3});
        val = str2double(tk{4});
        
        names{n} = sprintf('v%d%d%d', i, j, k_val);
        idx_i(n) = i;
        idx_j(n) = j;
        idx_k(n) = k_val;
        values(n) = val;
    end
    
    % 创建 Table 并排序
    T = table(names, idx_i, idx_j, idx_k, values, ...
        'VariableNames', {'Name', 'i', 'j', 'k', 'Value'});
    
    % 按 i (盐度), j (温度), k (压力) 的顺序排序
    T_sorted = sortrows(T, {'i', 'j', 'k'});
    
    % 4. 格式化打印
    fprintf('%-8s | %3s %3s %3s | %-22s\n', 'Name', 'i', 'j', 'k', 'Coefficient Value');
    fprintf('---------|-------------|-----------------------\n');
    
    for row = 1:height(T_sorted)
        fprintf('%-8s | %3d %3d %3d | %22.15e\n', ...
            T_sorted.Name{row}, ...
            T_sorted.i(row), ...
            T_sorted.j(row), ...
            T_sorted.k(row), ...
            T_sorted.Value(row));
    end
    
    fprintf('---------|-------------|-----------------------\n');
    fprintf('共计提取并打印: %d 项系数\n', height(T_sorted));
    
    if height(T_sorted) == 75
        fprintf('[检查通过] 系数数量正确 (75项)。\n');
    else
        fprintf('[检查警告] 系数数量异常 (预期75项)。\n');
    end
end