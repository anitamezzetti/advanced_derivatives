% Advanced Derivatives - problem set 5
% Hien Le, Francesco Maizza, Anita Mezzetti

clc
clear all

global lenk

% Inizialization: 

data = xlsread('SX5E_Impliedvols.xlsx');  %load data

S0 = 2772.70;               % current spot given by the paper

K = data(:,1) * S0;         % strikes from 40% to 200%
K = K(~isnan(K));           % delete NaN values

T = data(1,:);              % Maturity dates
T = T(~isnan(T));           % delete NaN values

lenk = length(K);           % number of strikes
lenT = length(T);           % number of maturity dates

dK = K(2)-K(1);             % vector dk(i) = k(i+1)-k(i)  
dT = [T(1),diff(T)];        % vector dT(i) = T(i+1)-T(i)  

% sempifications:
r = 0;                      % interst rate
q = 0;                      % dividend

vol = data(2:end,2:end);    % volatility data
positive_vol = vol>0;       % true/false matrix: positive => 1, zero => 0

vol_tilde = diag(K)*vol;    % vol_tilde(K,T) = K * vol(K,T) 

% Call Prices: algorithm to extrapolate call prices using the BS formula:
call_obs = zeros(lenk, lenT); % initialise the observed call prices matrix

for i = 1:lenk
    for j = 1:lenT
        if (positive_vol(i,j) >0) % otherwhise is not a stochastic process
            % find the call price for maturity T(j) and strike K(i):
            call_obs(i,j) = CallBS(S0,K(i),r,T(j),vol(i,j),q);
        end
    end
end

% Implied Volatility for the expirations in the spreadsheet:

C0 = max(S0-K,0);                   % initial call price for each K
result_C = zeros(lenk, lenT);        % resulting C after the AH algorithm
esimated_vol = zeros(size(vol));    % estimated volatility through interp. 

for j = 1:length(T)
    
    display(j)
    
    C_actual = C0;                              
    
    % find nonzero vol for all vol tilde of that maturity (one row for each T):
    pos_vol_T = find((vol_tilde(:,j)>0));
    
    %first guess corresponding to the observed sigma volatilities
    para0 = nonzeros(vol_tilde(:,j)); 
    
    %display(length(para0))
    
    % sigma lower and upper bounds for the optimization method
    lb = zeros(1,length(para0));        % sigma > 0
    up = S0 * ones(1,length(para0));    % sigma < S0
    
    %dt = dT(j);
    %dk = dK(j);
    
    
    % Optimization
    [parameters,fval,~,~] = fminsearchcon (@(parameters)optimization_function...
        (C_actual,call_obs(:,j),dT(j),dK,parameters),para0,lb,up);
    
    % Prediction (similar steps to optimization_function)
    pos_call_obs = find(call_obs(:,j)); % positions in which call_obs not zero
    
    % estimated volatility
    interp_vol = volatility_interpolation (pos_call_obs, parameters);
    esimated_vol(:,j) = interp_vol;
    
    inv_A = pinv(build_A(interp_vol,dT(j),dK)); % A^(-1) = psudoinverse of A
    % inv_A>=0 (showed by Nabben) imples that the discrete system is stable
    if inv_A<0
        fprintf("The discrete system is not stable")
    end
    C_next = inv_A * C_actual;   
    C_actual = C_next;      % update the currenr call price
    result_C(:,j) = C_next; % update result_C
    
end

% Prices for T=1 and T=1.5
TT = [1, 1.5];
C_models = zeros(lenk,length(TT)); % resulting C

for i = 1:length(TT)
    
    % first time T exceeds TT(i) minus 1 => last time T lower that TT(i)
    j = find(T>TT(i),1)-1;       %T_prime between T_j and T_{j+1}
    dT = TT(i)-T(j);
    sigma = esimated_vol(:,j);   % volatilities of that moment 
    inv_A = pinv(build_A(sigma,dT,dK));
    C_models(:,i) = inv_A * result_C(:,j);
end

% Unique computed call prices matrix:
% C_total is the resulting matrix of all call prices, at the beginning is
% composed by only result_C, then we add the last two columns for TT
% T_toal contains T and TT in the right order

% C_total = result_C;
% T_total = T;
% for i=1:length(TT)
%     
%     % j is needed to understand in what point we add the column for TT(i)
%     % call prices
%     j = find(T>TT(i),1)-1;
%     T_total = [T_total(1:j), TT(i), T_total(j+1:end)]; %add TT(i)
%     C_total=[C_models_total(:,1:j),C_models(:,i),C_models_total(:,j+1:end)];
% end


% Implied Volatility: (from the computed call prices)
sigma0 = 0.02;                  % initial value, pr 0.05
[IV,fval2] = implied_volBS(S0,K,r,T,q,result_C,sigma0);

% Volatility Surface Plot:

%kk = repmat(K,1,lenT);
%tt = repmat(T(2:end),lenk,1);

figure
surf(T,K,IV)
title('Volatility surface')
xlabel('Strikes')
ylabel('Maturities')
zlabel('Volatility')

% hold on
% % Observed values:
% [X,Y]=find(vol>0); 
% Z=vol(vol>0); 
% plot3(T(Y),K(X),Z,'.r','markersize',10)
% hold off

% Functions:

function f = optimization_function (C_actual,call_obs,dT,dK,para)
% return the function to minimize
% C0 = initial call price
% C_actual = option price ad current step
% call_obs = observed call prices for that T
% lenk = number of strikes and length of C0

    
pos_call_obs = find(call_obs); % positions in which call_obs not zero

% Volatility Interpolation:
interpolated_vol = volatility_interpolation (pos_call_obs, para);

% C(T(j+1)) = A^(-1) * C(T(j))
inv_A = pinv(build_A(interpolated_vol,dT,dK)); % A^(-1) = psudoinverse of A

% inv_A>=0 (showed by Nabben) imples that the discrete system is stable
if inv_A<0
    fprintf("The discrete system is not stable")
end
%display(inv_A)
%display(C_actual)
%display(size(inv_A))
%display(size(C_actual))

C_next = inv_A * C_actual;              % call price next step

% function to minimize: 
f=sum((C_next(pos_call_obs) - call_obs(pos_call_obs,1)).^2);
end

function sigma = volatility_interpolation (positions, sigma_in)
% returns interpolated values of a 1-D function at specific query points 
% using linear interpolation 
% positions: sample pointsc
% sigma_est: corresponding values of the sample points
% sigma: result

global lenk

% Interpolation of sigmas between the volatilities. 
sigma = interp1(positions, sigma_in, (1: lenk)', 'nearest');
% 'nearest': nearest neighbor interpolation
% For a point where the volatility is now known yet, it uses the 
% nearest observed value as image, which is in line with the problem.

% the length of sigma is always equal to the number of strikes (lenk)

% all the sigmas before the first estimated sigma set to the first
% estimated
sigma(1:positions(1)) = sigma(positions(1));
% all the sigmas after the last estimated sigma set to the last
% estimated
sigma(positions(end):end) = sigma(positions(end));

if length(sigma) ~= lenk
    fprintf("Error in the interpolation part");
end

end

function A = build_A (sigma,dT,dK)
% A is a tri-diagonal matrix, diagonally dominante with positive diagonal
% and negative off-diagonals

% define Z as written in explained in slide 12 of lecture 4
Z = 0.5 * dT / (dK)^2 * sigma(2:end).^2;% definition of Z
Z = Z(1:end-1)';                            % skip last value


diag_0 = [1, 1+2*Z, 1];                     % principal diagonal
diag_up1 = [0, -Z];                         % upper diagonal
diag_less1 = [-Z, 0];                       % lower digonal

ones_diag_0 = ones(length(diag_0),1);       % vector of ones for the diag

A = diag(diag_0, 0) + diag(diag_up1,1) + diag(diag_less1,-1);
end




