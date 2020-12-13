%% parameters
% for swaption
t = 0; % first date
r = 0.05; % term structure (flat)
T_s = 1; % first ex date
T_x = 2; % last ex date

% for underlying swap
tau = 0.25;
T_N = 4; % last payment date
payment_dates = (0:16).*tau;

% corresponding forward rates
vols_fw = zeros(16,1);
vols_fw(2:7) = 0.2;
vols_fw(8:11) = .22;
vols_fw(12:16) = .24;

thetas = pi/2 * ((1:16) - 2)./14;

% other params
N_sim = 10000;
K = .05;
N_x = length(tau:tau:T_x);
N_n = length(tau:tau:T_N);
Hs = 1e-4:1e-4:r;

%% computation for best H

[F,S] = ForwardRates(r,vols_fw,thetas,K,tau,N_sim,N_x,N_n);
%disp(S)
first_ex_date = T_s/tau;
S = S(:,first_ex_date+1:end);
%disp(S)
[~,num_ex_dates] = size(S); 
best_Hs = zeros(num_ex_dates-1,1);

% price of swaption at 4 dates
for t=1:4
    avg_swaption_prices = zeros(length(Hs),1);
    S_t = S(:,end-t:end);
    
    % discounted factor (slide 12 Lecture 11)
    D_t = cumprod(1./(1+tau*F(:,end-t:end,end-t)),2);
    %disp(D_t);
    %if (t==4)
    %    disp(D_t)
    %end
    
    % experimenting with different H
    for i=1:length(Hs)
        [~,dim2_S] = size(S_t);
        H = Hs(i).*ones(dim2_S,1);
        %disp(H)
        CF_H = Swaptions(S_t,H,D_t);
        %disp(size(D_t))
        avg_swaption_prices(i) = mean(CF_H); 
    end
    [~,argmax] = max(avg_swaption_prices);
    best_Hs(t) = Hs(argmax);
end

%% Generate new paths and compute prices with optimal Hs
[F,S] = ForwardRates(r,vols_fw,thetas,K,tau,N_sim,N_x,N_n);

first_ex_date = T_s/tau;
S = S(:,first_ex_date:end);
[~,num_ex_dates] = size(S); 

SPs = S;
S_t = S(:,end-4:end);
D_t = cumprod(1./(1+tau*F(:,end-4:end,end-t)),2);
CF_t = Swaptions(S_t,best_Hs,D_t);
fprintf("The swaption price is %.7f", mean(CF_t))
%The swaption price is 0.0053037


% 
% for t=1:4
%     S_t = S(:,end-t:end);
%     
%     % discounted factor
%     D_t = cumprod(1./(1+tau*F(:,end-t:end,end-t)),2);
% 
%     CF_t = Swaptions(S_t,best_Hs(t),D_t);
%     %disp(mean(CF_t));
%     SPs(:,t) = (CF_t);
% end

%disp("Price of Swaption is", mean(SPs)*10000);
