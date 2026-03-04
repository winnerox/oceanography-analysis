%% Check_EN4_Units.m
clear; clc;

Files = {'EN4_TSLA_Terms_1to15_Safe.mat', 'EN4_TSLA_Terms_1to15_GPU.mat', 'EN4_TSLA_Terms_1to15_StdRef.mat'};
TargetFile = '';
for i = 1:length(Files)
    if exist(Files{i}, 'file'), TargetFile = Files{i}; break; end
end

if isempty(TargetFile)
    fprintf('No EN4 file found.\n');
    return;
end

fprintf('Loading %s ...\n', TargetFile);
load(TargetFile, 'TSLA_AllOrders');

% Check Term 1 values (Thermal Expansion)
% Expected: 10-100 mm (if mm), or 0.01-0.1 (if m)
T1 = TSLA_AllOrders(:,:,:,1);
MeanVal = nanmean(abs(T1(:)));
MaxVal = nanmax(abs(T1(:)));

fprintf('Term 1 Stats:\n');
fprintf('  Mean Abs: %.4f\n', MeanVal);
fprintf('  Max Abs:  %.4f\n', MaxVal);

if MeanVal < 1
    fprintf('>> CONCLUSION: Unit is Likely METERS.\n');
else
    fprintf('>> CONCLUSION: Unit is Likely MILLIMETERS.\n');
end
