% 检查 Ishii CrossDetail 文件结构
fprintf('检查 Ishii_CrossDetail_Avg.mat 文件结构...\n');
load('d:\work\Ishii_mat_data\Ishii_CrossDetail_Avg.mat');

% 显示变量列表
fprintf('\n文件中的变量：\n');
vars = fieldnames(who);
for i = 1:length(vars)
    fprintf('%s\n', vars{i});
end

% 显示关键变量的大小
fprintf('\n关键变量大小：\n');
if exist('lon', 'var'), fprintf('lon: %s\n', mat2str(size(lon))); end
if exist('lat', 'var'), fprintf('lat: %s\n', mat2str(size(lat))); end
if exist('time_vec', 'var'), fprintf('time_vec: %s\n', mat2str(size(time_vec))); end
if exist('Cross_T1S1', 'var'), fprintf('Cross_T1S1: %s\n', mat2str(size(Cross_T1S1))); end

fprintf('\n检查完成！\n');