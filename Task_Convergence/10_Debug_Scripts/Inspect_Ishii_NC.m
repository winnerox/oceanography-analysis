%% Inspect_Ishii_NC.m
File = 'D:\work\Ishii_05_24\Temperature\temp.2005.nc';
vi = ncinfo(File);
disp('Variables:');
for i=1:length(vi.Variables)
    v=vi.Variables(i);
    fprintf('  %s: Size=[%s]\n', v.Name, num2str(v.Size));
end
