%--------------------------------------------------------------------------
% Advanced Derivatives 
% Problem Set 5
% Marc-Aurèle Divernois
%
% This code uses SX5E_Impliedvols.xlsx and AndreasenHuge.m and
% BlackScholes.m
%--------------------------------------------------------------------------
% Exercice 1
%--------------------------------------------------------------------------
clear all;
clc;

% Extract data

data = xlsread('SX5E_Impliedvols.xlsx'); % col1 : SPX, col2 : AMZN
T = [0 data(1,2:end)];
ImpVols = data(2:end,2:end);
t=0;
r=0;
q=0;
S0=2772.70;
K = data(2:end,1)*S0;
type = "call";

% Get option prices from Implied Vols

LiGrid = length(K);
CoGrid = length(T);
C = NaN(LiGrid,CoGrid);
if(type == "call")
    C0 = max(S0-K,0)';
elseif(type == "put")
    C0 = max(K-S0,0)';
else
    print('type error')
end

C(:,1) = C0;

for i=1:LiGrid
    for j =1:CoGrid-1
        C(i,j+1) = BlackScholes(S0,t,ImpVols(i,j),K(i),T(j+1),r,q,type);
    end
end

IndexesToConsider = ImpVols>0;
C(:,2:end) = C(:,2:end).*IndexesToConsider;

% Optimization problem

for j = 1:CoGrid-1 
    NoNaNindex = find((ImpVols(:,j)>0));
    Nparam = length(NoNaNindex);
    deltaT = T(j+1)-T(j);
    deltaK = K(2)-K(1);
    fun = @(v) sum(IndexesToConsider(:,j).*(AndreasenHuge(NoNaNindex,j,v,C,deltaT,deltaK,K)-C(:,j+1)).^2);
    init = 0.2*ones(length(NoNaNindex),1)*S0;
    LB = zeros(Nparam,1); 
    UB = S0*ones(Nparam,1);
    [V{j},fval] = fmincon(fun,init,[],[],[],[],LB,UB);
    %[V{j},fval] = fminsearch(fun,init);

    InterpolatedC = AndreasenHuge(NoNaNindex,j,V{j},C,deltaT,deltaK,K);
    C(:,j+1) = InterpolatedC;
end

% Get implied vols from interpolated prices

Interpol_IV = NaN(size(ImpVols));

for i = 1:LiGrid
    for j = 2:CoGrid
        fun2 = @(sigma) (BlackScholes(S0,t,sigma,K(i),T(j),r,q,type) - C(i,j))^2;
        [Interpol_IV(i,j-1),fval(i,j)] = fminsearch(fun2,0.52);
    end
end


% Price surface

x = repmat(K,1,CoGrid);
y = repmat(T,LiGrid,1);

figure
surf(x,y,C)
title('Call prices surface')
xlabel('Strikes')
ylabel('Maturities')
zlabel('Price')


% Volatility surface

x = repmat(K,1,CoGrid-1);
y = repmat(T(2:end),LiGrid,1);

figure
surf(x,y,Interpol_IV)
title('Volatility surface')
xlabel('Strikes')
ylabel('Maturities')
zlabel('Volatility')

% Extrapolation
% Expiration T=1 and T=1.5

wantedT = 1;
%wantedT = 1.5;
j = find(T>1,1)-1;
NoNaNindex = find((ImpVols(:,j)>0));
Nparam = length(NoNaNindex);
deltaT = wantedT-T(j);
deltaK = K(2)-K(1);

AH = AndreasenHuge(NoNaNindex,j,V{j-1},C,deltaT,deltaK,K);




