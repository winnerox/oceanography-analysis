% 1. 定义文件夹路径 (请确保这是解压后的文件夹)
folderPath = 'D:\科研\EN4_analyses_c13_last20years\EN.4.2.2.analyses.c13.2005';

% 2. 搜索文件夹下的所有 .nc 文件
fileList = dir(fullfile(folderPath, '*.nc'));

% 3. 检查是否找到了文件
if isempty(fileList)
    error('在该文件夹内未找到 .nc 文件。请确认：1. 路径是否正确；2. 压缩包是否已解压。');
else
    % 获取第一个文件的完整路径
    targetFile = fullfile(fileList(1).folder, fileList(1).name);
    fprintf('正在读取文件: %s\n', targetFile);
    
    % 使用 ncdisp 查看变量
    ncdisp(targetFile);
end