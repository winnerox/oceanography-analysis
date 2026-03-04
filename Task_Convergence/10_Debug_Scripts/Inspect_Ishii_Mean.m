%% Inspect_Ishii_Mean.m
clear; clc;
LoadFile = 'D:\work\Ishii_TSLA_Terms\Ishii_Mean_State.mat';

if ~exist(LoadFile, 'file')
    error('File not found: %s', LoadFile);
end

S = load(LoadFile);
Vars = fieldnames(S);

fprintf('File: %s\n', LoadFile);
fprintf('Variables:\n');
for i = 1:length(Vars)
    val = S.(Vars{i});
    fprintf('  %s: %s, Class: %s\n', Vars{i}, mat2str(size(val)), class(val));
    
    if isnumeric(val) && numel(val) > 1
        fprintf('     Min: %f, Max: %f, NaNs: %d\n', min(val(:)), max(val(:)), sum(isnan(val(:))));
        if all(val(:) == 0)
            fprintf('     [WARNING] All Zeros!\n');
        end
    end
end
