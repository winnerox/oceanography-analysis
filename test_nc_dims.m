clear; clc;
try
    disp('--- IAP TEMP ---');
    info = ncinfo('D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc');
    for i = 1:length(info.Variables)
        if strcmp(info.Variables(i).Name, 'temp') || strcmp(info.Variables(i).Name, 'salinity')
            fprintf('%s: ', info.Variables(i).Name);
            disp(info.Variables(i).Size);
        end
    end
    
    disp('--- IAP SALT ---');
    info = ncinfo('D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc');
    for i = 1:length(info.Variables)
        if strcmp(info.Variables(i).Name, 'temp') || strcmp(info.Variables(i).Name, 'salinity')
            fprintf('%s: ', info.Variables(i).Name);
            disp(info.Variables(i).Size);
        end
    end
    
    disp('--- Ishii SALT ---');
    info = ncinfo('D:\work\Ishii_05_24\Salinity\sal.2005.nc');
    for i = 1:length(info.Variables)
        if strcmp(info.Variables(i).Name, 'sal') || strcmp(info.Variables(i).Name, 'salinity')
            fprintf('%s: ', info.Variables(i).Name);
            disp(info.Variables(i).Size);
        end
    end
catch ME
    disp(ME.message);
end
exit;
