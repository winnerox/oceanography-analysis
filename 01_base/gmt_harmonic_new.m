function [ Amplitude1, Amplitude1_std, Phase1, Phase1_std, Amplitude2, Amplitude2_std, Phase2, Phase2_std, Trend, Trend_std, Trend_sig, Trend_p, Trend_t, Trend_line, Resid, Interp] = gmt_harmonic_new(t,t1,grid_data,grid_data_std)
%
% Copyright: FENG Wei, Institute of Geodesy and Geophysics, Chinese Academy
% of Sciences
%
% FENG Wei 25/03/2015

% INPUT
% t              is time epoch
% t1             is the time epoch for the predicted (missing) months
% grid_data      is the input field grid_data(lat,lon,epochs), 2-D or 3-D matrix
% grid_data_std  is the standard deviations of input field

% OUTPUT
% Amplitude1     = amplitude of annual cycle
% Amplitude1_std = standard deviation of annual amplitude
% Phase1         = phase of annual cycle, angle of cosine
% Phase1_std     = standard deviation of annual phase
% Amplitude2     = amplitude of semi-annual cycle
% Amplitude2_std = standard deviation of semi-annual amplitude
% Phase2         = phase of semi-annual cycle
% Phase2_std     = standard deviation of semi-annual phase
% Trend          = trend of time series
% Trend_std      = standard deviation of estimated trend
% Trend_sig      = trend significance (1=significant, 0=not significant)
% Trend_p        = p-value of trend
% Trend_t        = t-statistic of trend
% Trend_line     = predicted trend line
% Interp         = predicted trend + seasonal cycles
% Resid          = residuals after removing trend and seasonal cycles

if ndims(grid_data)==2  % grid is 2-D matrix
    rows=1;
    cols=1;
    epochs= max(size(t));                             %     epochs= max(size(grid_data)); change by zyl
    grid_data=reshape(grid_data,[rows,cols,epochs]);
    if nargin<3 || nargin>4
        error('Number of input parameters in gmt_harmonic function is wrong!');
    end
    if nargin==3
        grid_data_std = [];
    end
end

[rows,cols,epochs] = size(grid_data);

if nargin<3 || nargin>4
    error('Number of input parameters in gmt_harmonic function is wrong!');
end

if nargin==3
    grid_data_std = [];
else
    grid_data_std=reshape(grid_data_std,[rows,cols,epochs]);
end

omega = 2*pi;

% ------------------
% initial
% ------------------
time_series     = zeros(1,epochs);
Trend           = zeros(rows, cols);
Trend_std       = zeros(rows, cols);
Trend_sig       = zeros(rows, cols); % 趋势显著性 (1=显著, 0=不显著)
Trend_p         = zeros(rows, cols); % 趋势p值
Trend_t         = zeros(rows, cols); % 趋势t统计量
Amplitude1      = zeros(rows, cols);
Amplitude1_std  = zeros(rows, cols);
Phase1          = zeros(rows, cols);
Phase1_std      = zeros(rows, cols);
Amplitude2      = zeros(rows, cols);
Amplitude2_std  = zeros(rows, cols);
Phase2          = zeros(rows, cols);
Phase2_std      = zeros(rows, cols);
Trend_line      = zeros(rows, cols, epochs);
Resid           = zeros(rows, cols, epochs);
Interp          = zeros(rows, cols, epochs);
% ------------------
% end
% ------------------


% create normal matrix
A = [];
% ones function create unit matrix with epochs X 1
A(:,1) = (ones(epochs,1));
A(:,2) = t';
A(:,3) = cos(omega*t'); % annual
A(:,4) = sin(omega*t');
A(:,5) = cos(omega*2*t'); % semi-annual
A(:,6) = sin(omega*2*t');
NM = (A'*A)\A';%最小二乘的表达式

% No standard deviations of input field
if nargin==3
    for i = 1:rows
        for j = 1:cols
            for k = 1:epochs
                time_series(k) = grid_data(i,j,k);
            end
            if sum(isnan(time_series))==0 % there is no NaN in time series
                % normal equation, x contains all estimated parameters
                x    = NM*time_series';
                % cofactor matrix
                Exx = (A'*A)\eye(size(A'*A));
                %Exx = inv((A'*A));
                
                %
                % annual ampltidue & phase
                %
                Ampl = sqrt(x(3)^2 + x(4)^2);
                if x(4)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(3)/Ampl);
                elseif x(4)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(3)/Ampl);
                elseif x(4)==0 && x(3)==1
                    Pha = 0;
                elseif x(4)==0 && x(3)==-1
                    Pha = pi;
                elseif x(4)==0 && x(3)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360;
                Amplitude1(i,j) = Ampl; % annual amplitude
                Phase1(i,j)     = Pha; % annual phase
                
                %
                % semi-annual ampltidue & phase
                %
                Ampl = sqrt(x(5)^2 + x(6)^2);
                if x(6)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(5)/Ampl);
                elseif x(6)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(5)/Ampl);
                elseif x(6)==0 && x(5)==1
                    Pha = 0;
                elseif x(6)==0 && x(5)==-1
                    Pha = pi;
                elseif x(6)==0 && x(5)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360/2;
                Amplitude2(i,j) = Ampl; % semi-annual amplitude
                Phase2(i,j)      = Pha; % semi-annual phase
                
                % predict time series at epochs t1 using bias, trend, seasonal
                % cycles
                if ~isempty(t1)
                    Trend_line(i,j,:) = x(1) + x(2)*t1';
                    Oscilation_int=  x(1) + x(2)*t1' + x(3)*cos(omega*t1') + x(4)*sin(omega*t1') + x(5)*cos(2*omega*t1') + x(6)*sin(2*omega*t1');
                    Interp(i,j,:)=Oscilation_int;
                end
                
                % residual
                res = time_series' - A*x;
                Resid(i,j,:) = res;
                
                % aposteriori variance estimate: Error^2/(number of obs - number of param)
                % aposteriori unit weight mean error 验后单位权中误差
                var_est = res'*res/(length(time_series) - 6);
                % aposteriori covariance matrix for estimated parameters from aposteriori
                % unit weight mean error and co-factor matrix
                % 协因数阵乘以单位权中误差得到验后的观测值权中误差矩阵
                Exx_scal = var_est*Exx;
                
                Trend(i,j)     = x(2);
                Trend_std(i,j) = 2*sqrt(Exx_scal(2,2));
                
                % ========== 新增：趋势显著性检验 ==========
                % 计算t统计量
                Trend_t(i,j) = x(2) / sqrt(Exx_scal(2,2));
                
                % 计算p值 (双尾检验)
                degrees_freedom = epochs - 6; % 自由度 = 观测数 - 参数个数
                Trend_p(i,j) = 2 * (1 - tcdf(abs(Trend_t(i,j)), degrees_freedom));
                
                % 判断显著性 (5%显著性水平)
                Trend_sig(i,j) = Trend_p(i,j) < 0.05;
                % ========== 结束新增 ==========
                
                % 1 sigma: 68.3% ; 2 sigma: 95.4%; 3 sigma: 99.7%
                Amplitude1_std(i,j) = 2*sqrt( ( x(3)*x(3) * Exx_scal(3,3) + x(4)*x(4) * Exx_scal(4,4) )/(Amplitude1(i,j)*Amplitude1(i,j)));
                if abs(x(4)/x(3))<=1
                    x_43 = x(4)/x(3);
                    Phase1_std(i,j) = (1-x_43^2+x_43^4-x_43^6+x_43^8-x_43^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_34 = x(3)/x(4);
                    Phase1_std(i,j) =( x_34^2-x_34^4+x_34^6-x_34^8+x_34^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
                
                Amplitude2_std(i,j) = 2*sqrt( ( x(5)*x(5) * Exx_scal(5,5) + x(6) * x(6) * Exx_scal(6,6) )/(Amplitude2(i,j)*Amplitude2(i,j)));
                if abs(x(6)/x(5))<=1
                    x_65 = x(6)/x(5);
                    Phase2_std(i,j) = (1-x_65^2+x_65^4-x_65^6+x_65^8-x_65^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_56 = x(5)/x(6);
                    Phase2_std(i,j) =( x_56^2-x_56^4+x_56^6-x_56^8+x_56^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
            else % there is NaN in time series
                Amplitude1(i,j) = NaN; % annual amplitude
                Amplitude1_std(i,j)=NaN;
                Phase1(i,j)     = NaN; % annual phase
                Phase1_std(i,j) = NaN;
                Amplitude2(i,j) = NaN; % semi-annual amplitude
                Phase2(i,j)      = NaN; % semi-annual phase
                Phase2_std(i,j) = NaN;
                Trend(i,j)     = NaN;
                Trend_std(i,j) = NaN;
                Trend_sig(i,j) = NaN; % 新增
                Trend_p(i,j)   = NaN; % 新增
                Trend_t(i,j)   = NaN; % 新增
                Trend_line(i,j,1:epochs)=NaN;
                Resid(i,j,1:epochs) =NaN;
                Interp(i,j,1:epochs)=NaN;
            end
        end
    end
    
    % There are standard deviations of input field
elseif nargin==4
    for i = 1:rows
        for j = 1:cols
            for k = 1:epochs
                time_series(k) = grid_data(i,j,k);
                time_series_var(k) = grid_data_std(i,j,k)^2;
            end
            if sum(isnan(time_series))==0 % there is no NaN in time series
                % set a priori unit weight mean error is ONE
                % weights
                P = diag(1./time_series_var);
                
                % normal equation
                % x    = (A'*P*A)\eye(size(A'*P*A))*A'*P*time_series';
                x    = (A'*P*A)\A'*P*time_series';
                
                %  cofactor matrix
                %Exx = (A'*P*A)\eye(size(A'*P*A));
                Exx = inv((A'*P*A));
                
                %
                % annual ampltidue & phase
                %
                Ampl = sqrt(x(3)^2 + x(4)^2);
                if x(4)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(3)/Ampl);
                elseif x(4)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(3)/Ampl);
                elseif x(4)==0 && x(3)==1
                    Pha = 0;
                elseif x(4)==0 && x(3)==-1
                    Pha = pi;
                elseif x(4)==0 && x(3)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360;
                Amplitude1(i,j) = Ampl;
                Phase1(i,j)      = Pha;
                
                %
                % semi-annual  ampltidue & phase
                %
                Ampl = sqrt(x(5)^2 + x(6)^2);
                if x(6)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(5)/Ampl);
                elseif x(6)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(5)/Ampl);
                elseif x(6)==0 && x(5)==1
                    Pha = 0;
                elseif x(6)==0 && x(5)==-1
                    Pha = pi;
                elseif x(6)==0 && x(5)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360/2;
                Amplitude2(i,j) = Ampl;
                Phase2(i,j)      = Pha;
                
                
                % predict time series at epochs t1 using bias, trend, seasonal
                % cycles
                if ~isempty(t1)
                    Trend_line(i,j,:) = x(1) + x(2)*t1';
                    Oscilation_int=  x(1) + x(2)*t1' + x(3)*cos(omega*t1') + x(4)*sin(omega*t1') + x(5)*cos(2*omega*t1') + x(6)*sin(2*omega*t1');
                    Interp(i,j,:)=Oscilation_int;
                end
                
                % residual
                res = time_series' - A*x;
                Resid(i,j,:) = res;
                
                %2020-8-28做了修改，误差肯定是验后的，所以。。。。。
%                 % covariance matrix for estimated parameters
%                 Exx_scal = Exx; % 用验前的单位权中误差            


                % aposteriori variance estimate: Error^2/(number of obs - number of param)
                % aposteriori unit weight mean error 验后单位权中误差
                var_est = res'*res/(length(time_series) - 6);
                % aposteriori covariance matrix for estimated parameters from aposteriori
                % unit weight mean error and co-factor matrix
                % 协因数阵乘以单位权中误差得到验后的观测值权中误差矩阵
                Exx_scal = var_est*Exx;
                
                Trend(i,j)     = x(2);
                Trend_std(i,j) = 2*sqrt(Exx_scal(2,2));
                
                % ========== 新增：趋势显著性检验 ==========
                % 计算t统计量
                Trend_t(i,j) = x(2) / sqrt(Exx_scal(2,2));
                
                % 计算p值 (双尾检验)
                degrees_freedom = epochs - 6; % 自由度 = 观测数 - 参数个数
                Trend_p(i,j) = 2 * (1 - tcdf(abs(Trend_t(i,j)), degrees_freedom));
                
                % 判断显著性 (5%显著性水平)
                Trend_sig(i,j) = Trend_p(i,j) < 0.05;
                % ========== 结束新增 ==========
                   
                % 1 sigma: 68.3% ; 2 sigma: 95.4%; 3 sigma: 99.7%
                Amplitude1_std(i,j) = 2*sqrt( ( x(3)*x(3) * Exx_scal(3,3) + x(4)*x(4) * Exx_scal(4,4) )/(Amplitude1(i,j)*Amplitude1(i,j)));
                if abs(x(4)/x(3))<=1
                    x_43 = x(4)/x(3);
                    Phase1_std(i,j) = (1-x_43^2+x_43^4-x_43^6+x_43^8-x_43^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_34 = x(3)/x(4);
                    Phase1_std(i,j) =( x_34^2-x_34^4+x_34^6-x_34^8+x_34^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
                
                Amplitude2_std(i,j) = 2*sqrt( ( x(5)*x(5) * Exx_scal(5,5) + x(6) * x(6) * Exx_scal(6,6) )/(Amplitude2(i,j)*Amplitude2(i,j)));
                if abs(x(6)/x(5))<=1
                    x_65 = x(6)/x(5);
                    Phase2_std(i,j) = (1-x_65^2+x_65^4-x_65^6+x_65^8-x_65^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_56 = x(5)/x(6);
                    Phase2_std(i,j) =( x_56^2-x_56^4+x_56^6-x_56^8+x_56^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
            else % there is NaN in time series
                Amplitude1(i,j) = NaN; % annual amplitude
                Amplitude1_std(i,j)=NaN;
                Phase1(i,j)     = NaN; % annual phase
                Phase1_std(i,j) = NaN;
                Amplitude2(i,j) = NaN; % semi-annual amplitude
                Phase2(i,j)      = NaN; % semi-annual phase
                Phase2_std(i,j) = NaN;
                Trend(i,j)     = NaN;
                Trend_std(i,j) = NaN;
                Trend_sig(i,j) = NaN; % 新增
                Trend_p(i,j)   = NaN; % 新增
                Trend_t(i,j)   = NaN; % 新增
                Trend_line(i,j,1:epochs)=NaN;
                Resid(i,j,1:epochs) =NaN;
                Interp(i,j,1:epochs)=NaN;
            end
        end
    end
end