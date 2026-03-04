%% Ishii 数据处理完整方案
clear; clc;
addpath('D:\work');

%% 1. 基础配置
DataRoot = 'D:\work\Ishii_05_24';
TempDir = fullfile(DataRoot, 'Temperature');
SaltDir = fullfile(DataRoot, 'Salinity');
OutputDir = 'D:\work\Ishii_TSLA_Terms';

Years = 2005:2024;
Months = 1:12;
MaxOrder = 8;
MaxDepth = 2000;
rho0 = single(1035.0);

if ~exist(OutputDir, '