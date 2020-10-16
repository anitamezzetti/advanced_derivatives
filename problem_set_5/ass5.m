% Advanced Derivatives - problem set 5
% Hien Le, Francesco Maizza, Anita Mezzetti

clc
clear all
close all

% Inizialization: 

data = xlsread('SX5E_Impliedvols.xlsx');  %load data

S0 = 2772.70; %bohhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

K = data(:,1) * S0;         % strikes from 40% to 200%
K = K(~isnan(K));           % delete NaN values

T = data(1,:);              % Maturity dates
T = T(~isnan(T));           % delete NaN values

lenk = length(K);          % numbe of strikes
lenT = length(T);          % number of maturity dates

dk = K(2:end) - K(1:end-1); % vector dk(i) = k(i+1)-k(i)  
dT = T(2:end) - T(1:end-1); % vector dT(i) = T(i+1)-T(i)  

% semifications:
r = 0;                      % interst rate
q = 0;                      % dividend

vol = data(2:end,2:end);
vol_marked = 1.000000000000;

% Call Prices: algorithm to extrapolate call prices using the BS formula:
call_obs = zeros(lenk, lenT); % initialise the observed call prices matrix

for i = 1:lenk
    for j = 1:lenT
        % find the call price for maturity T(j) and strike K(i)
        call_obs(i,j) = CallBS(S0,K(i),r,T(j),vol(i,j),q);
        
        % not sureeeeeeeeeeeeeeeee 
        if call_obs(i,j)<0
            call_obs(i,j) = 0;

        end
    end
end

% Implied Volatility for the expirations in the spreadsheet:
for j = 1:lenT
    
    
    
    
end

% Prices for T=1 and T=1.5
TT = [1, 1.5];




% Volatility Surface




% Plot


