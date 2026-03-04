clc
clear
close all

%海水比容、热容和盐容计算程序

%这是目前的最新版本=========================
%%这里的平均态包含三类==============
%1)05-15年，多年平均
%2)标准的温度和盐度35，0
%3)05-24年，多年平均
Datasets='EN4_c13';



%需要手动修改的
start_year      = 2005;
start_mon       = 1;
end_year        = 2024;
end_mon         = 12;
Dep=2000;  %SSLA,TSLA,HSLA 计算的深度 0-Dep的结果


%初始化参数  
rho_0 = 1028.0;  %海水密度

num = (end_year-start_year)*12 + end_mon - (start_mon-1); % number of months

start_time = start_year+start_mon/12-1/24;
end_time = end_year+end_mon/12-1/24;
time1 = start_time:1/12:end_time;

int_year=floor(time1);
int_mon=round(((time1-int_year)*24+1)/2);
time_argo = int_year+int_mon./12-1/24;
time_dim = length(time_argo);

%读取基本信息========================================================
file = 'H:\ARGO\EN4_c13\EN.4.2.2.f.analysis.c13.200501.nc';

depth1 = ncread(file,'depth');   %积分层的平均深度
lat=(-90:89)';
lon=(1:360)';
lon_dim = length(lon);
lat_dim = length(lat);
dep_dim = length(depth1);
T_mean = NaN(lon_dim,lat_dim,dep_dim);
S_mean = NaN(lon_dim,lat_dim,dep_dim);
depth_bnds = ncread(file,'depth_bnds');
depth2 = depth_bnds(2,:)';   %积分层的底部深度
        
density = NaN(lon_dim,lat_dim,dep_dim);    
beta = NaN(lon_dim,lat_dim,dep_dim);    
alpha = NaN(lon_dim,lat_dim,dep_dim);

T_all = NaN(lon_dim,lat_dim,dep_dim,num);
S_all = NaN(lon_dim,lat_dim,dep_dim,num);
T = NaN(lon_dim,lat_dim,dep_dim);       %温度异常 
S = NaN(lon_dim,lat_dim,dep_dim);       %盐度异常

% %计算depth 和 dz=========================================================
lat_grid = repmat(lat',lon_dim,1,dep_dim);
lon_grid = repmat(lon,1,lat_dim,dep_dim);

depth = repmat(depth1,1,lon_dim,lat_dim);
depth = permute(depth,[2 3 1]);
press = gsw_p_from_z(-depth ,lat_grid);  %深度换算成压强gsw

%计算不同的积分dz区间=======================================================
tmp = depth_bnds(2,:)-depth_bnds(1,:);
dz_all = repmat(tmp',1,lon_dim,lat_dim);
dz_all = permute(dz_all,[2 3 1]);

%==========================================================================
cumsum_dz = cumsum(dz_all,3);
cumsum_dz(cumsum_dz> 2000)= 2000;
cumsum_dz(:,:,2:end)=cumsum_dz(:,:,2:end)-cumsum_dz(:,:,1:end-1);%先累加，后替换，再累减，省去了插值的过程。
dz_2000a = cumsum_dz;    %2000m以上       
dz_2000b = dz_all -dz_2000a; %2000m以下


%计算平均态密度===================================================


for kk=1:num
   
	dir_in_full = strcat('H:\ARGO\EN4_c13\EN.4.2.2.f.analysis.c13.',num2str(int_year(kk)),num2str(int_mon(kk),'%02u'),'.nc');%赋值dir_in_full，文件名赋值
	if exist(dir_in_full,'file')            %dir_in_full只有41个字符 char， 

        T_all(:,8:180,:,kk) = ncread(dir_in_full,'temperature')-273.15;
        S_all(:,8:180,:,kk) = ncread(dir_in_full,'salinity');
    end
    
end


%提取0m——2000m——4000m的掩膜====================
T_tmp1 = T_all(:,:,:,1);
T_Nan = isnan(T_tmp1);
S_tmp1 = T_all(:,:,:,1);
S_Nan = isnan(S_tmp1);

T_tmp = T_all(:,:,1,1);
S_tmp = S_all(:,:,1,1);
index_0_new=isnan(T_tmp);
index_2000 = find(depth2 >= 2000, 1, 'first');
T_tmp = T_all(:,:,index_2000,1);
T_nan_2000_new = permute(T_tmp,[2 1]);
index_2000_new=isnan(T_nan_2000_new);

index_4000 = find(depth2 >= 4000, 1, 'first');
T_tmp = T_all(:,:,index_4000,1);
T_nan_4000_new = permute(T_tmp,[2 1]);
index_4000_new=isnan(T_nan_4000_new);


%计算平均态密度=======================================================
index_05_15 = find(time_argo>2005 & time_argo<2016); 
index_05_24 = find(time_argo>2005 & time_argo<2025);

S_mean_05_15 = mean(S_all(:,:,:,index_05_15), 4);
T_mean_05_15 = mean(T_all(:,:,:,index_05_15), 4);
S_mean_05_15 = gsw_SA_from_SP(S_mean_05_15,press,lon_grid,lat_grid);       %SP转为SA
T_mean_05_15 = gsw_CT_from_pt(S_mean_05_15,T_mean_05_15);    %位温转化为CT
D_mean_05_15 = gsw_rho_CT_exact(S_mean_05_15,T_mean_05_15,press); %计算平均态下的密度 gsw

S_mean_05_24 = mean(S_all(:,:,:,index_05_24), 4);
T_mean_05_24 = mean(T_all(:,:,:,index_05_24), 4);
S_mean_05_24 = gsw_SA_from_SP(S_mean_05_24,press,lon_grid,lat_grid);       %SP转为SA
T_mean_05_24 = gsw_CT_from_pt(S_mean_05_24,T_mean_05_24);    %位温转化为CT
D_mean_05_24 = gsw_rho_CT_exact(S_mean_05_24,T_mean_05_24,press); %计算平均态下的密度 gsw

%计算参考态密度=======================================================
S_ref = 35.16504.* ones(lon_dim,lat_dim,dep_dim);
T_ref = zeros(lon_dim,lat_dim,dep_dim);

T_ref(T_Nan) = NaN;
S_ref(S_Nan) = NaN;

%计算比容、热容和盐容====================================================
SSLA = NaN(lon_dim,lat_dim,time_dim,15);
TSLA = NaN(lon_dim,lat_dim,time_dim,15);
HSLA = NaN(lon_dim,lat_dim,time_dim,15);
SSLA_error = NaN(lon_dim,lat_dim,time_dim,9);

SSLA_2000a = NaN(lon_dim,lat_dim,time_dim,15);
TSLA_2000a = NaN(lon_dim,lat_dim,time_dim,15);
HSLA_2000a = NaN(lon_dim,lat_dim,time_dim,15);
SSLA_2000a_error = NaN(lon_dim,time_dim,dep_dim,9);

SSLA_2000b = NaN(lon_dim,lat_dim,time_dim,15);
TSLA_2000b = NaN(lon_dim,lat_dim,time_dim,15);
HSLA_2000b = NaN(lon_dim,lat_dim,time_dim,15);
SSLA_2000b_error = NaN(lon_dim,lat_dim,time_dim,9);

for kk=1:num
    tic
	dir_in_full = strcat('H:\ARGO\EN4_c13\EN.4.2.2.f.analysis.c13.',num2str(int_year(kk)),num2str(int_mon(kk),'%02u'),'.nc');%赋值dir_in_full，文件名赋值
	if exist(dir_in_full,'file')            %dir_in_full只有41个字符 char， 

        T(:,8:180,:) = ncread(dir_in_full,'temperature')-273.15;
        S(:,8:180,:)  = ncread(dir_in_full,'salinity');
        S = gsw_SA_from_SP(S, press,lon_grid,lat_grid);       %SP转为SA
        T = gsw_CT_from_pt(S,T);    %位温转化为CT

        rho = gsw_rho_CT_exact(S,T,press); %gsw
        grav = gsw_grav(lat_grid,press); %计算重力加速度
        rho_a = rho(:,:,1);
        grav_a = grav(:,:,1);
        rho_m = nanmean(rho,3);
        %====================================================================================================
        %计算比容=========平均密度为准05-15
        rho_ref = D_mean_05_15; %gsw 

        SSLA(:,:,kk,1) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        SSLA(:,:,kk,2) = -nansum((rho-rho_ref).*dz_all,3)./rho_0;
        SSLA(:,:,kk,3) = -nansum((rho-rho_ref).*dz_all,3)./rho_m;
        SSLA(:,:,kk,4) = -nansum((rho-rho_ref).*dz_all,3)./rho_a;
        SSLA(:,:,kk,5) = -nansum((rho-rho_ref).*dz_all.*grav,3)./(rho_a.*grav_a);
        
        %计算比容=========参考密度为准
        rho_ref = gsw_rho_CT_exact(S_ref,T_ref,press); %gsw 

        SSLA(:,:,kk,6) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        SSLA(:,:,kk,7) = -nansum((rho-rho_ref).*dz_all,3)./rho_0;
        SSLA(:,:,kk,8) = -nansum((rho-rho_ref).*dz_all,3)./rho_m;
        SSLA(:,:,kk,9) = -nansum((rho-rho_ref).*dz_all,3)./rho_a;
        SSLA(:,:,kk,10) = -nansum((rho-rho_ref).*dz_all.*grav,3)./(rho_a.*grav_a);

        %计算比容=========平均密度为准05-24
        rho_ref = D_mean_05_24; %gsw 
        
        SSLA(:,:,kk,11) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        SSLA(:,:,kk,12) = -nansum((rho-rho_ref).*dz_all,3)./rho_0;
        SSLA(:,:,kk,13) = -nansum((rho-rho_ref).*dz_all,3)./rho_m;
        SSLA(:,:,kk,14) = -nansum((rho-rho_ref).*dz_all,3)./rho_a;
        SSLA(:,:,kk,15) = -nansum((rho-rho_ref).*dz_all.*grav,3)./(rho_a.*grav_a);

        %=========================================================================
        %计算热容和盐容=========平均密度为准05-15=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T_mean_05_15,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T_mean_05_15,press);
        T_anomoly = T - T_mean_05_15;
        S_anomoly = S - S_mean_05_15;

        TSLA(:,:,kk,1) = nansum(alpha.*T_anomoly.*dz_all,3);
        HSLA(:,:,kk,1) = -nansum(beta.*S_anomoly.*dz_all,3);

        TSLA(:,:,kk,2) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        HSLA(:,:,kk,2) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,1) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_all,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA(:,:,kk,4) = nansum(alpha.*T_anomoly.*dz_all,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA(:,:,kk,4) = -nansum(beta.*S_anomoly.*dz_all,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_15,press);
        TSLA(:,:,kk,5) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,2) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T,press);
        HSLA(:,:,kk,5) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,3) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);
        %====================================================================================================    

        %计算热容和盐容=========参考密度为准=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T_ref,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T_ref,press);
        T_anomoly = T - T_ref;
        S_anomoly = S - S_ref;

        TSLA(:,:,kk,6) = nansum(alpha.*T_anomoly.*dz_all,3);
        HSLA(:,:,kk,6) = -nansum(beta.*S_anomoly.*dz_all,3);

        TSLA(:,:,kk,7) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        HSLA(:,:,kk,7) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,4) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_ref,press); %gsw
        TSLA(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        rho_ref = gsw_rho_CT_exact(S_ref,T,press); %gsw
        HSLA(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_all,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_ref,press); %gsw
        TSLA(:,:,kk,9) = nansum(alpha.*T_anomoly.*dz_all,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T,press); %gsw
        HSLA(:,:,kk,9) = -nansum(beta.*S_anomoly.*dz_all,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_ref,press);
        TSLA(:,:,kk,10) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,5) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T,press);
        HSLA(:,:,kk,10) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,6) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);

        
        %计算热容和盐容=========平均密度为准05-24=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T_mean_05_24,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T_mean_05_24,press);
        T_anomoly = T - T_mean_05_24;
        S_anomoly = S - S_mean_05_24;

        TSLA(:,:,kk,11) = nansum(alpha.*T_anomoly.*dz_all,3);
        HSLA(:,:,kk,11) = -nansum(beta.*S_anomoly.*dz_all,3);

        TSLA(:,:,kk,12) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        HSLA(:,:,kk,12) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,7) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);

        rho_ref = gsw_rho_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_all,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_all,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA(:,:,kk,14) = nansum(alpha.*T_anomoly.*dz_all,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA(:,:,kk,14) = -nansum(beta.*S_anomoly.*dz_all,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_24,press);
        TSLA(:,:,kk,15) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,8) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T,press);
        HSLA(:,:,kk,15) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_all,3);
        SSLA_error(:,:,kk,9) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_all,3);
        %==================================================================================================== 
        %     2000a
        %     2000a
        %====================================================================================================
        %====================================================================================================
        %计算比容=========平均密度为准05-15
        rho_ref = D_mean_05_15; %gsw 

        SSLA_2000a(:,:,kk,1) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        SSLA_2000a(:,:,kk,2) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_0;
        SSLA_2000a(:,:,kk,3) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_m;
        SSLA_2000a(:,:,kk,4) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_a;
        SSLA_2000a(:,:,kk,5) = -nansum((rho-rho_ref).*dz_2000a.*grav,3)./(rho_a.*grav_a);
        
        %计算比容=========参考密度为准
        rho_ref = gsw_rho_CT_exact(S_ref,T_ref,press); %gsw 

        SSLA_2000a(:,:,kk,6) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        SSLA_2000a(:,:,kk,7) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_0;
        SSLA_2000a(:,:,kk,8) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_m;
        SSLA_2000a(:,:,kk,9) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_a;
        SSLA_2000a(:,:,kk,10) = -nansum((rho-rho_ref).*dz_2000a.*grav,3)./(rho_a.*grav_a);

        %计算比容=========平均密度为准05-24
        rho_ref = D_mean_05_24; %gsw 
        
        SSLA_2000a(:,:,kk,11) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        SSLA_2000a(:,:,kk,12) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_0;
        SSLA_2000a(:,:,kk,13) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_m;
        SSLA_2000a(:,:,kk,14) = -nansum((rho-rho_ref).*dz_2000a,3)./rho_a;
        SSLA_2000a(:,:,kk,15) = -nansum((rho-rho_ref).*dz_2000a.*grav,3)./(rho_a.*grav_a);

        %=========================================================================
        %计算热容和盐容=========平均密度为准05-15=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T_mean_05_15,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T_mean_05_15,press);
        T_anomoly = T - T_mean_05_15;
        S_anomoly = S - S_mean_05_15;

        TSLA_2000a(:,:,kk,1) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        HSLA_2000a(:,:,kk,1) = -nansum(beta.*S_anomoly.*dz_2000a,3);

        TSLA_2000a(:,:,kk,2) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        HSLA_2000a(:,:,kk,2) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,1) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA_2000a(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA_2000a(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA_2000a(:,:,kk,4) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA_2000a(:,:,kk,4) = -nansum(beta.*S_anomoly.*dz_2000a,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_15,press);
        TSLA_2000a(:,:,kk,5) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,2) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T,press);
        HSLA_2000a(:,:,kk,5) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,3) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);
        %====================================================================================================    

        %计算热容和盐容=========参考密度为准=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T_ref,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T_ref,press);
        T_anomoly = T - T_ref;
        S_anomoly = S - S_ref;

        TSLA_2000a(:,:,kk,6) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        HSLA_2000a(:,:,kk,6) = -nansum(beta.*S_anomoly.*dz_2000a,3);

        TSLA_2000a(:,:,kk,7) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        HSLA_2000a(:,:,kk,7) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,4) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_ref,press); %gsw
        TSLA_2000a(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        rho_ref = gsw_rho_CT_exact(S_ref,T,press); %gsw
        HSLA_2000a(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_ref,press); %gsw
        TSLA_2000a(:,:,kk,9) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T,press); %gsw
        HSLA_2000a(:,:,kk,9) = -nansum(beta.*S_anomoly.*dz_2000a,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_ref,press);
        TSLA_2000a(:,:,kk,10) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,5) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T,press);
        HSLA_2000a(:,:,kk,10) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,6) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);

        
        %计算热容和盐容=========平均密度为准05-24=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T_mean_05_24,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T_mean_05_24,press);
        T_anomoly = T - T_mean_05_24;
        S_anomoly = S - S_mean_05_24;

        TSLA_2000a(:,:,kk,11) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        HSLA_2000a(:,:,kk,11) = -nansum(beta.*S_anomoly.*dz_2000a,3);

        TSLA_2000a(:,:,kk,12) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        HSLA_2000a(:,:,kk,12) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,7) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);

        rho_ref = gsw_rho_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA_2000a(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA_2000a(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_2000a,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA_2000a(:,:,kk,14) = nansum(alpha.*T_anomoly.*dz_2000a,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA_2000a(:,:,kk,14) = -nansum(beta.*S_anomoly.*dz_2000a,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_24,press);
        TSLA_2000a(:,:,kk,15) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,8) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T,press);
        HSLA_2000a(:,:,kk,15) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000a,3);
        SSLA_2000a_error(:,:,kk,9) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000a,3);
        %====================================================================================================   
        %==================================================================================================== 
        %     2000b
        %     2000b
        %====================================================================================================
        %====================================================================================================
        %====================================================================================================
        %计算比容=========平均密度为准05-15r
        rho_ref = D_mean_05_15; %gsw 

        SSLA_2000b(:,:,kk,1) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        SSLA_2000b(:,:,kk,2) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_0;
        SSLA_2000b(:,:,kk,3) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_m;
        SSLA_2000b(:,:,kk,4) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_a;
        SSLA_2000b(:,:,kk,5) = -nansum((rho-rho_ref).*dz_2000b.*grav,3)./(rho_a.*grav_a);
        
        %计算比容=========参考密度为准
        rho_ref = gsw_rho_CT_exact(S_ref,T_ref,press); %gsw 

        SSLA_2000b(:,:,kk,6) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        SSLA_2000b(:,:,kk,7) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_0;
        SSLA_2000b(:,:,kk,8) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_m;
        SSLA_2000b(:,:,kk,9) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_a;
        SSLA_2000b(:,:,kk,10) = -nansum((rho-rho_ref).*dz_2000b.*grav,3)./(rho_a.*grav_a);

        %计算比容=========平均密度为准05-24
        rho_ref = D_mean_05_24; %gsw 
        
        SSLA_2000b(:,:,kk,11) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        SSLA_2000b(:,:,kk,12) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_0;
        SSLA_2000b(:,:,kk,13) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_m;
        SSLA_2000b(:,:,kk,14) = -nansum((rho-rho_ref).*dz_2000b,3)./rho_a;
        SSLA_2000b(:,:,kk,15) = -nansum((rho-rho_ref).*dz_2000b.*grav,3)./(rho_a.*grav_a);

        %=========================================================================
        %计算热容和盐容=========平均密度为准05-15=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T_mean_05_15,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T_mean_05_15,press);
        T_anomoly = T - T_mean_05_15;
        S_anomoly = S - S_mean_05_15;

        TSLA_2000b(:,:,kk,1) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        HSLA_2000b(:,:,kk,1) = -nansum(beta.*S_anomoly.*dz_2000b,3);

        TSLA_2000b(:,:,kk,2) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        HSLA_2000b(:,:,kk,2) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,1) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA_2000b(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA_2000b(:,:,kk,3) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_15,press); %gsw
        TSLA_2000b(:,:,kk,4) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_15,T,press); %gsw
        HSLA_2000b(:,:,kk,4) = -nansum(beta.*S_anomoly.*dz_2000b,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_15,press);
        TSLA_2000b(:,:,kk,5) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,2) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_15,T,press);
        HSLA_2000b(:,:,kk,5) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,3) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);
        %====================================================================================================    

        %计算热容和盐容=========参考密度为准=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T_ref,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T_ref,press);
        T_anomoly = T - T_ref;
        S_anomoly = S - S_ref;

        TSLA_2000b(:,:,kk,6) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        HSLA_2000b(:,:,kk,6) = -nansum(beta.*S_anomoly.*dz_2000b,3);

        TSLA_2000b(:,:,kk,7) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        HSLA_2000b(:,:,kk,7) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,4) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);

        
        rho_ref = gsw_rho_CT_exact(S,T_ref,press); %gsw
        TSLA_2000b(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        rho_ref = gsw_rho_CT_exact(S_ref,T,press); %gsw
        HSLA_2000b(:,:,kk,8) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_ref,press); %gsw
        TSLA_2000b(:,:,kk,9) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_ref,T,press); %gsw
        HSLA_2000b(:,:,kk,9) = -nansum(beta.*S_anomoly.*dz_2000b,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_ref,press);
        TSLA_2000b(:,:,kk,10) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,5) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_ref,T,press);
        HSLA_2000b(:,:,kk,10) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,6) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);

        
        %计算热容和盐容=========平均密度为准05-24=====================
        [rho_ref, alpha, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T_mean_05_24,press); %gsw
        [rho_SA_SA, rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T_mean_05_24,press);
        T_anomoly = T - T_mean_05_24;
        S_anomoly = S - S_mean_05_24;

        TSLA_2000b(:,:,kk,11) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        HSLA_2000b(:,:,kk,11) = -nansum(beta.*S_anomoly.*dz_2000b,3);

        TSLA_2000b(:,:,kk,12) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        HSLA_2000b(:,:,kk,12) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,7) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);

        rho_ref = gsw_rho_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA_2000b(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 
        rho_ref = gsw_rho_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA_2000b(:,:,kk,13) = nansum((1.0-rho./rho_ref).*dz_2000b,3); 

        [~, alpha, ~] = gsw_rho_alpha_beta_CT_exact(S,T_mean_05_24,press); %gsw
        TSLA_2000b(:,:,kk,14) = nansum(alpha.*T_anomoly.*dz_2000b,3);
        [~, ~, beta] = gsw_rho_alpha_beta_CT_exact(S_mean_05_24,T,press); %gsw
        HSLA_2000b(:,:,kk,14) = -nansum(beta.*S_anomoly.*dz_2000b,3); 

        [~,rho_SA_CT, rho_CT_CT, ~, ~] = gsw_rho_second_derivatives_CT_exact(S,T_mean_05_24,press);
        TSLA_2000b(:,:,kk,15) = nansum((alpha.*(T_anomoly) - 0.5.*T_anomoly.*T_anomoly.*rho_CT_CT./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,8) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);
        [rho_SA_SA, rho_SA_CT, ~, ~, ~] = gsw_rho_second_derivatives_CT_exact(S_mean_05_24,T,press);
        HSLA_2000b(:,:,kk,15) = -nansum((beta.*S_anomoly + 0.5.*S_anomoly.*S_anomoly.*rho_SA_SA./rho_ref).*dz_2000b,3);
        SSLA_2000b_error(:,:,kk,9) = -nansum((rho_SA_CT.* S_anomoly .* T_anomoly./rho_ref).*dz_2000b,3);

    end
    toc
    disp([num2str(kk),'/',num2str(num)]);
    
end
%计算比容、热容和盐容====================================================
save([Datasets '_all.mat'],'SSLA','TSLA','HSLA','SSLA_error','SSLA_2000a','TSLA_2000a','HSLA_2000a','SSLA_2000a_error',...
    'SSLA_2000b','TSLA_2000b','HSLA_2000b','SSLA_2000b_error','index_2000_new','index_4000_new','time_argo','lat','lon');


