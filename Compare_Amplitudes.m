%% Compare_Amplitudes.m
clear; clc;

% --- 1. EN4 ---
files = {'EN4_TSLA_Terms_1to15_Safe.mat', 'EN4_TSLA_Terms_1to15_GPU.mat', 'EN4_TSLA_Terms_1to15_StdRef.mat'};
Target = '';
for i=1:length(files), if exist(files{i}, 'file'), Target=files{i}; break; end; end
if ~isempty(Target)
    fprintf('Loading EN4: %s\n', Target);
    load(Target, 'TSLA_AllOrders');
    T1 = squeeze(TSLA_AllOrders(:,:,:,1));
    % Smart Unit Detection for EN4
    if nanmean(abs(T1(:))) < 5, F=1000; else, F=1; end
    
    % Calc Amp
    Tmp = T1(1:10:end); % Subsample for speed
    Tmp = Tmp(~isnan(Tmp));
    % Valid points only logic approx
    % Just take std of whole vector vs time? No, need point-wise std.
    % Do simplified check: Global STD of T1 anomaly
    % Better: Calculate Amp for first 100 valid points
    [r,c] = find(~isnan(TSLA_AllOrders(:,:,1,1)));
    Vals = [];
    CT = 0;
    for k=1:min(length(r), 100)
       ts = squeeze(TSLA_AllOrders(r(k),c(k),:,1));
       if sum(~isnan(ts)) > 100
           p = polyfit(1:length(ts), ts', 1);
           detrend = ts - polyval(p, (1:length(ts))');
           Vals(end+1) = std(detrend) * F;
       end
    end
    fprintf('EN4  Median T1 Amp: %.4f mm\n', median(Vals));
else
    fprintf('EN4 Not Found\n');
end

% --- 2. IAP ---
try
    load('D:\work\MAT_Data\IAP_TSLA_Terms_1to8_Average.mat', 'TSLA_AllOrders');
    T1 = squeeze(TSLA_AllOrders(:,:,:,1));
     if nanmean(abs(T1(:))) < 5, F=1000; else, F=1; end
    
    [r,c] = find(~isnan(TSLA_AllOrders(:,:,1,1)));
    Vals = [];
    for k=1:min(length(r), 100)
       ts = squeeze(TSLA_AllOrders(r(k),c(k),:,1));
       if sum(~isnan(ts)) > 100
           p = polyfit(1:length(ts), ts', 1);
           detrend = ts - polyval(p, (1:length(ts))');
           Vals(end+1) = std(detrend) * F;
       end
    end
    fprintf('IAP  Median T1 Amp: %.4f mm\n', median(Vals));
catch
    fprintf('IAP Not Found\n');
end

% --- 3. Ishii ---
try
    load('D:\work\MAT_Data\Ishii_TSLA_Terms_1to8_Average.mat', 'TSLA_AllOrders');
    T1 = squeeze(TSLA_AllOrders(:,:,:,1));
     if nanmean(abs(T1(:))) < 5, F=1000; else, F=1; end
    
    [r,c] = find(~isnan(TSLA_AllOrders(:,:,1,1)));
    Vals = [];
    for k=1:min(length(r), 100)
       ts = squeeze(TSLA_AllOrders(r(k),c(k),:,1));
       if sum(~isnan(ts)) > 100
           p = polyfit(1:length(ts), ts', 1);
           detrend = ts - polyval(p, (1:length(ts))');
           Vals(end+1) = std(detrend) * F;
       end
    end
    fprintf('Ishii Median T1 Amp: %.4f mm\n', median(Vals));
catch
    fprintf('Ishii Not Found\n');
end
