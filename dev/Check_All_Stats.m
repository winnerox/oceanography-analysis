%% Check_All_Stats.m
clear; clc;

% List files
Files = {
    'EN4',   'EN4_TSLA_Terms_1to15_StdRef.mat';
    'IAP',   'D:\work\MAT_Data\IAP_TSLA_Terms_1to8_Average.mat';
    'Ishii', 'D:\work\MAT_Data\Ishii_TSLA_Terms_1to8_Average.mat';
};

fprintf('%-10s | %-15s | %-15s | %-15s\n', 'Model', 'Mean(Abs)', 'STD(Global)', 'UnitEst');
fprintf('%s\n', repmat('-',1,60));

for i = 1:size(Files,1)
    Name = Files{i,1};
    Path = Files{i,2};
    
    if ~exist(Path, 'file')
        % Try finding EN4 in current dir
        if strcmp(Name, 'EN4')
            d = dir('EN4_TSLA_Terms_1to15_*.mat');
            if ~isempty(d), Path = d(1).name; end
        end
    end
    
    if exist(Path, 'file')
        try
            dat = load(Path, 'TSLA_AllOrders');
            T1 = squeeze(dat.TSLA_AllOrders(:,:,:,1));
            
            Mu = nanmean(abs(T1(:)));
            Sigma = nanstd(T1(:));
            
            if Mu < 5, Unit='m'; Scale=1000; else, Unit='mm'; Scale=1; end
            
            fprintf('%-10s | %10.4f (%s) | %10.4f (%s) | %10.4f (mm)\n', ...
                Name, Mu, Unit, Sigma, Unit, Sigma * Scale);
        catch
             fprintf('%-10s | Error Loading\n', Name);
        end
    else
        fprintf('%-10s | Not Found\n', Name);
    end
end
