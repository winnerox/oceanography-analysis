%% Calc_CrossDetail_Ishii_Avg.m
% =========================================================================
% 功能：计算 Ishii 的详细混合项 (Order 2-8) [平均态]
% 输出：D:\work\Ishii_TSLA_Terms\Ishii_CrossDetail_Avg.mat
% =========================================================================
clear; clc;

DataDir = 'D:\work\Ishii_05_24';
OutputDir = 'D:\work\Ishii_TSLA_Terms';
MeanFile = fullfile(OutputDir, 'Ishii_Mean_State.mat');
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end

Years = 2005:2024; rho0 = 1035.0; MaxDepth = 2000; MaxOrder = 8;
CacheFile = fullfile('D:\work\EN4_TSLA_Terms', 'TEOS10_Engine_Cache.mat');
Engine = TEOS10_HighOrder_Engine(MaxOrder, CacheFile, false);

if ~exist(MeanFile, 'file'), error('请先运行 Calc_Mean_Ishii.m'); end
load(MeanFile, 'Mean_T', 'Mean_S', 'Lon', 'Lat', 'Depth', 'DepthIdx');

CalcDepth = Depth(DepthIdx);

% Fix Lon/Lat for GSW - Robust Cleaning
Lon = double(Lon); Lat = double(Lat);
if any(isnan(Lon)), Lon(isnan(Lon)) = 0; end
if any(isnan(Lat)), Lat(isnan(Lat)) = 0; end
Lon = mod(Lon, 360); 
Lat(Lat < -90) = -90; Lat(Lat > 90) = 90;
Depth = double(Depth); % Force double

[LON_3D, LAT_3D, DEPTH_3D] = ndgrid(Lon, Lat, CalcDepth);
P_3D = gsw_p_from_z(-DEPTH_3D, LAT_3D);
[Nx, Ny, Nz] = size(P_3D);

dz = zeros(Nz, 1);
dz(1) = 0.5 * (CalcDepth(1) + CalcDepth(2));
for k = 2:Nz-1, dz(k) = 0.5 * (CalcDepth(k+1) - CalcDepth(k-1)); end
dz(Nz) = CalcDepth(Nz) - CalcDepth(Nz-1);
dz_perm = reshape(dz, 1, 1, []);

% Derivatives
fprintf('Calculating Derivatives...\n');
% Robust SA_Mean Calculation
SA_Mean = zeros(size(P_3D));
for k = 1:Nz
    try
        SA_Mean(:,:,k) = gsw_SA_from_SP(squeeze(Mean_S(:,:,k)), squeeze(P_3D(:,:,k)), ...
                                        squeeze(LON_3D(:,:,k)), squeeze(LAT_3D(:,:,k)));
    catch
        SA_Mean(:,:,k) = squeeze(Mean_S(:,:,k)) * (35.16504 / 35.0);
    end
end
CT_Mean = gsw_CT_from_pt(SA_Mean, Mean_T);
Derivs = cell(MaxOrder, MaxOrder+1);

for n = 2:MaxOrder
    for k = 1:n-1
        n_T = n - k;
        D_val = Engine.calculate_mixed(SA_Mean, CT_Mean, P_3D, n_T, k);
        Derivs{n, k+1} = single(D_val);
    end
end

Variables = {};
for n = 2:MaxOrder
    for k = 1:n-1
       VarName = sprintf('Cross_T%dS%d', n-k, k);
       Variables{end+1} = VarName;
       eval([VarName ' = zeros(Nx, Ny, length(Years)*12, ''single'');']);
    end
end

InvFact = zeros(21,1); for i=0:20, InvFact(i+1)=1/factorial(i); end
TotalSteps = length(Years) * 12;
time_vec = zeros(TotalSteps, 1);
idx = 0;

for y = Years
    try
        PatT = sprintf('temp.%d.nc', y);
        PatS = sprintf('sal.%d.nc', y);
        FT = fullfile(DataDir, 'Temperature', PatT);
        FS = fullfile(DataDir, 'Salinity', PatS);
        
        if ~exist(FT, 'file') || ~exist(FS, 'file'), continue; end
        
        % Check var names
        vi = ncinfo(FT); vnames = {vi.Variables.Name};
        if any(strcmp(vnames, 'temp')), vT='temp'; elseif any(strcmp(vnames, 'temperature')), vT='temperature'; elseif any(strcmp(vnames, 'pt')), vT='pt'; else, vT=vnames{end}; end
        vi = ncinfo(FS); vnames = {vi.Variables.Name};
        if any(strcmp(vnames, 'salt')), vS='salt'; elseif any(strcmp(vnames, 'salinity')), vS='salinity'; elseif any(strcmp(vnames, 'so')), vS='so'; else, vS=vnames{end}; end

        for m = 1:12
            idx = idx + 1;
            time_vec(idx) = y + (m-0.5)/12;

            T_raw = double(ncread(FT, vT, [1 1 1 m], [Inf Inf Nz 1]));
            S_raw = double(ncread(FS, vS, [1 1 1 m], [Inf Inf Nz 1]));
            
            T_raw(T_raw > 100) = NaN;
            S_raw(S_raw > 100) = NaN;

            if nanmean(T_raw(:)) > 100, T_raw = T_raw - 273.15; end
            
            % Robust SA
            SA = zeros(size(S_raw));
            for k = 1:Nz
                 try
                     SA(:,:,k) = gsw_SA_from_SP(squeeze(S_raw(:,:,k)), squeeze(P_3D(:,:,k)), ...
                                            squeeze(LON_3D(:,:,k)), squeeze(LAT_3D(:,:,k)));
                 catch
                     SA(:,:,k) = squeeze(S_raw(:,:,k)) * (35.16504 / 35.0);
                 end
            end
            
            CT = gsw_CT_from_pt(SA, T_raw);
            
            dS = SA - SA_Mean; dT = CT - CT_Mean;
            
            S_Pow = ones(Nx, Ny, Nz, MaxOrder+1); T_Pow = ones(Nx, Ny, Nz, MaxOrder+1);
            for p = 1:MaxOrder
                S_Pow(:,:,:,p+1) = S_Pow(:,:,:,p) .* dS;
                T_Pow(:,:,:,p+1) = T_Pow(:,:,:,p) .* dT;
            end
            
            for n = 2:MaxOrder
                for k = 1:n-1
                    nT = n - k; nS = k;
                    C = InvFact(nS+1) * InvFact(nT+1);
                    term = C * Derivs{n, k+1} .* S_Pow(:,:,:,nS+1) .* T_Pow(:,:,:,nT+1);
                    val = -(nansum(term .* dz_perm, 3) / rho0) * 1000;
                    VarName = sprintf('Cross_T%dS%d', nT, nS);
                    eval([VarName '(:,:,idx) = single(val);']);
                end
            end
        end
        fprintf('Year %d Done.\n', y);
        
    catch ME
        fprintf('   [Error] Year %d: %s\n', y, ME.message);
    end
end

SaveName = fullfile(OutputDir, 'Ishii_CrossDetail_Avg.mat');
lon=Lon; lat=Lat;
SaveVars = {'lon', 'lat', 'time_vec'};
SaveVars = [SaveVars, Variables];
save(SaveName, SaveVars{:}, '-v7.3');
fprintf('Saved %s\n', SaveName);
