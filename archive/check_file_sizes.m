fprintf('=== 文件大小检查 ===\n\n');

folders = {'d:\work\EN4_TSLA_Terms', 'd:\work\IAP_TSLA_Terms', 'd:\work\Ishii_TSLA_Terms'};

for f = 1:length(folders)
    folder = folders{f};
    fprintf('【%s】\n', folder);
    
    files = dir(fullfile(folder, '*.mat'));
    fprintf('%-45s %10s\n', '文件名', '大小(MB)');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for i = 1:length(files)
        size_mb = files(i).bytes / (1024 * 1024);
        fprintf('%-45s %10.2f\n', files(i).name, size_mb);
    end
    fprintf('\n');
end
