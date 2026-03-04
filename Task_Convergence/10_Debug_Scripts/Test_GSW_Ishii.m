%% Test_GSW_Ishii.m
clear; clc;
DataDir = 'D:\work\Ishii_05_24';
List = dir(fullfile(DataDir, 'Temperature', '*.nc'));
SampleFile = fullfile(List(1).folder, List(1).name);

% 1. Load Dimensions (Same logic as failed script)
vi = ncinfo(SampleFile); vnames = {vi.Variables.Name};
if any(strcmp(vnames, 'lon')), vLon='lon'; elseif any(strcmp(vnames, 'long')), vLon='long'; elseif any(strcmp(vnames, 'longitude')), vLon='longitude'; else, vLon=vnames{1}; end
if any(strcmp(vnames, 'lat')), vLat='lat'; elseif any(strcmp(vnames, 'latitude')), vLat='latitude'; else, vLat=vnames{2}; end
if any(strcmp(vnames, 'depth')), vDep='depth'; elseif any(strcmp(vnames, 'depth_std')), vDep='depth_std'; else, vDep='z'; end

Lon = double(ncread(SampleFile, vLon));
Lat = double(ncread(SampleFile, vLat));
Depth = double(ncread(SampleFile, vDep));
MaxDepth = 2000;

% 2. Clean (Same as failed script)
Lon = mod(Lon, 360);
Lat(Lat < -89.9) = -89.9; Lat(Lat > 89.9) = 89.9;
if any(Depth < 0), Depth(Depth < 0) = 0; end

DepthIdx = find(Depth <= MaxDepth);
CalcDepth = Depth(DepthIdx);
[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
Std_S = ones(size(P_3D)) * 35.0;

fprintf('Dimensions: %d x %d x %d\n', length(Lon), length(Lat), length(CalcDepth));

% 3. Test Full Grid
try
    fprintf('Testing entire grid... ');
    SA = gsw_SA_from_SP(Std_S, P_3D, LON_3D, LAT_3D);
    fprintf('Success!\n');
catch ME
    fprintf('Failed: %s\n', ME.message);
    
    % 4. Isolate Error
    fprintf('Isolating error...\n');
    Nx = length(Lon); Ny = length(Lat);
    
    % Test Middle
    try
        gsw_SA_from_SP(35, P_3D(Nx/2, Ny/2, 1), LON_3D(Nx/2, Ny/2, 1), LAT_3D(Nx/2, Ny/2, 1));
        fprintf('  Single point (center) works.\n');
    catch, fprintf('  Single point (center) FAILED.\n'); end
    
    % Test -180..180 conversion
    try
        fprintf('Testing -180..180 conversion... ');
        Lon2 = Lon; Lon2(Lon2 > 180) = Lon2(Lon2 > 180) - 360;
        [L3, A3, ~] = ndgrid(Lon2, Lat, CalcDepth);
        SA = gsw_SA_from_SP(Std_S, P_3D, L3, A3);
        fprintf('Success!\n');
    catch ME2
        fprintf('Failed: %s\n', ME2.message);
    end
end
