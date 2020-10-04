clear all
clc

% set parameters
q = 0;
lambda = 1; %one jump per year
gamma = 0.1; % jump size
sigma = 0.2;
r = 0.04;
S0 = 100;

% expiration dates
T1 = 0.02;
T2 = 0.08;
T = [T1; T2];
t = 0;

K = linspace(60,150, 500)'; %strikes

C = zeros(length(K), 2); % initial price call

% discounted C with a jump process 
for i=1:length(K)       % for each strike
    for j=1:length(T)   % for each maturity date
        temp = 0;
        for k=1:100
            S0_bs = S0*(1 - gamma)^(k-1);
            q_bs = q-lambda*gamma;
            payoff_bs = black_scholes(S0_bs , t, sigma, K(i), T(j), r, q_bs);
            
            % probability k-1 jumps up to the maturity:
            jumps = exp(-lambda*T(j))*((lambda*T(j))^(k-1)/(factorial(k-1)));
            
            C(i,j) = C(i,j) + jumps * payoff_bs;
        end
    end
end

% undiscounted C
% Computing C tild the undiscouted price 
C_und = C.*exp(r*(T'-t));


% We compute phi by differenetation of C 
der_C = C_und(1:end-2,:) - 2 * C_und(2:end-1,:) + C_und(3:end,:);
phi = der_C/(K(2)-K(1))^2;

% We plot the result 
figure
plot(K(2:end-1,1), phi(:,1), 'g')
hold on
plot(K(2:end-1,1), phi(:,2), 'b')
hold off
xlabel('Strike')
ylabel('Phi')
title('Implied Probabilty Pistribution')
legend('maturity 0.02','maturity 0.08')


