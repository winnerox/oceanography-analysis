%% Inspect_IAP_NC.m
DataDir = 'D:\work\IAP_05_24\TEMP';
List = dir(fullfile(DataDir, '*.nc'));
if isempty(List)
    disp('No IAP files found.');
else
    File = fullfile(List(1).folder, List(1).name);
    disp(['Inspecting: ' File]);
    try
        info = ncinfo(File);
        disp('Variables:');
        for i = 1:length(info.Variables)
            disp(['  ' info.Variables(i).Name]);
        end
    catch ME
        disp(ME.message);
    end
end
