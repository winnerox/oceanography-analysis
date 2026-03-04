%% Check_Cross_Terms_Variables.m
% =========================================================================
% 功能：读取三个混合项文件，显示其中包含的变量
% =========================================================================
clear; clc;

% 定义输出文件
output_file = 'Cross_Terms_Variables.txt';
fid = fopen(output_file, 'w');

% 定义文件路径
files = {
    'd:\work\EN4_TSLA_Terms\EN4_Cross_Terms_1to8_Average.mat',
    'd:\work\IAP_mat_data\IAP_Cross_Terms_1to8_Average.mat',
    'd:\work\IAP_TSLA_Terms\IAP_Cross_Terms_1to8_Average.mat'
};

% 遍历每个文件
for i = 1:length(files)
    file_path = files{i};
    fprintf(fid, '\n=== 检查文件: %s ===\n', file_path);
    fprintf('\n=== 检查文件: %s ===\n', file_path);
    
    if exist(file_path, 'file')
        % 读取文件中的变量
        variables = whos('-file', file_path);
        
        if ~isempty(variables)
            fprintf(fid, '文件包含以下变量:\n');
            fprintf('文件包含以下变量:\n');
            for j = 1:length(variables)
                var = variables(j);
                fprintf(fid, '  - 名称: %s, 大小: %s, 类型: %s\n', ...
                    var.name, mat2str(var.size), var.class);
                fprintf('  - 名称: %s, 大小: %s, 类型: %s\n', ...
                    var.name, mat2str(var.size), var.class);
            end
        else
            fprintf(fid, '文件中没有变量\n');
            fprintf('文件中没有变量\n');
        end
    else
        fprintf(fid, '文件不存在\n');
        fprintf('文件不存在\n');
    end
end

fprintf(fid, '\n>>> 检查完成! <<<\n');
fprintf('\n>>> 检查完成! <<<\n');
fclose(fid);

fprintf('结果已保存到: %s\n', output_file);

