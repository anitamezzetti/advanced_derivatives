%--------------------------------------------------------------------------
% Advanced Derivatives 
% Problem Set 4
% Marc-Aurèle Divernois
%
% This code uses Impvols_SPX_AMZN.xlsx and BlackScholes.m
%--------------------------------------------------------------------------
% Exercice 1
%--------------------------------------------------------------------------
clear all;
clc;

% Extract data

data = xlsread('Impvols_SPX_AMZN.xlsx'); % col1 : SPX, col2 : AMZN
K = data(:,[1 5]);
ImpVol = data(:,[2 6]);
S0 = [2921 1971];
T = 108/365;
r = 0.024;
q = [1.8/100 1.9/100];
rho = 0.5;
t=0;

% Get option prices from implied vols

C = NaN(length(K),length(S0));
for i=1:length(K)
    for j =1:length(S0)
        C(i,j) = BlackScholes(S0(j),t,ImpVol(i,j),K(i,j),T,r,q(j));
    end
end

C = C*exp(r*(T-t)); % undiscounted price

% Implied CDF

Phi = NaN(length(K)-1,length(S0));
for j = 1:length(S0)
    Phi(:,j) = 1 + (C(2:end,j)-C(1:end-1,j))./(K(2:end,j)-K(1:end-1,j));
end
Phi(sum(~isnan(Phi(:,1))),1) = 1; %make sure that CDF ends at 1
Phi(sum(~isnan(Phi(:,2))),2) = 1;

% Generate correlated normally distributed random numbers

Nsim = 10000;
x = mvnrnd([0 0],[1 rho;rho 1],Nsim);
x1 = normcdf(x);

% Inverse transform sampling

index = NaN(size(x1));
S_simulated = NaN(Nsim,length(S0));
for i=1:Nsim
    for j=1:length(S0)
        index(i,j) = find(Phi(:,j)>x1(i,j),1);
        S_simulated(i,j) = K(index(i,j),j);
    end
end

% Compute payoffs and price of the outperformance option

payoff = exp(-r*(T-t))*max(0,S_simulated(:,1)/S0(1) - S_simulated(:,2)/S0(2));
price = mean(payoff);

% Price = [0.0555, 0.0585] depending on the simulation
