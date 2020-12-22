%clear all
%close all
%clc
analytical_price = 0.0074897997038153585;
%% params
K = 0.805; % strike
k = 0.15; % kappa
sigma_r = 0.01;
theta = 0.05;
r0 = 0.042;

T1 = 3/12; % 3 months
T2 = 5 + 3/12; % 5 years 3 months

%% other self-defined params
Nrs = 100:10:1000;
Nts = 5:5:150;
r_min = 0;
r_max = 1;

%% different rs
Vs = zeros(length(Nrs),1);
for i=1:length(Nrs)
    Nr = Nrs(i);
    %disp(Nrs(i))
    dt = 0.0001;
    Nt = T1/dt;
    Vs(i) = FiniteDiffCN(Nr, Nt, r_max, r_min, K, k, sigma_r, theta,...
        r0, T1, T2);
end

%% different ts
Vs_T = zeros(length(Nts),1); 
for i=1:length(Nts)
    Nt = Nts(i);
    %disp(Nrs(i))
    Nr = 500;
    Vs_T(i) = FiniteDiffCN(Nr, Nt, r_max, r_min, K, k, sigma_r, theta,...
        r0, T1, T2);
end

%% plot for Nr
figure(1)
plot(Nrs, Vs, 'r.-','LineWidth',2,'MarkerSize',10)
hold on
plot(Nrs, analytical_price*ones(length(Vs),1), 'b-', 'LineWidth',...
    2,'MarkerSize',10)
hold on
plot(Nrs, analytical_price*1.01*ones(length(Vs),1), 'k-', 'LineWidth',...
    0.75,'MarkerSize',5)
hold on
plot(Nrs, analytical_price*0.99*ones(length(Vs),1), 'y-', 'LineWidth',...
    0.75,'MarkerSize',5)
legend({'Crank-Nicholson', 'Analytical', 'Upper bound', 'Lower bound'})
xlabel("Nr")
ylabel("V0")

saveas(gcf,'CN_r.png')
hold off

%% plot for Nt
figure(2)
hold off
plot(Nts, Vs_T, 'r.-','LineWidth',2,'MarkerSize',10)
hold on
plot(Nts, 0.00748*ones(length(Vs_T),1), 'b-',...
    'LineWidth',2,'MarkerSize',10)
hold on
plot(Nts, analytical_price*1.01*ones(length(Vs_T),1), 'k-', 'LineWidth',...
    0.75,'MarkerSize',5)
hold on
plot(Nts, analytical_price*0.99*ones(length(Vs_T),1), 'y-', 'LineWidth',...
    0.75,'MarkerSize',5)

legend({'Crank-Nicholson', 'Analytical', 'Upper bound', 'Lower bound'})
xlabel("Nt")
ylabel("V0")

saveas(gcf,'CN_t.png')

%% Intervals where price is 1% from analytical
intervals_r = Nrs(abs());
intervals_t = Nts(abs(Vs_T - analytical_price) <= 0.01);

%% CN
function price = FiniteDiffCN(Nr, Nt, r_max, r_min, K, k, sigma_r,...
        theta, r0, T1, T2)

dr = (r_max-r_min)/Nr;
r_grid = r_min:dr:r_max; 

dt = T1/Nt;
t_grid = 0:dt:T1;

M = round((r_max - r_min)/dr); % number of rates

N= length(t_grid); % number of time points

V = zeros(M+1, N+1);
%disp(size(V))
% boundary conditions
V(:, N+1) = max(0, K - bond_price(T1, T2, sigma_r, k, theta, r_grid))'; 
V(:,1) = (K- bond_price(T1, T2, sigma_r, k, theta, 0))';
V(M+1, :) = 0;

% see overleaf file for explanation
a = @(r) (k*(theta-r))/(4*dr) + (sigma_r^2)/(4*dr^2);
b = @(r) (k*(theta-r))/(4*dr) - (sigma_r^2)/(4*dr^2);
c = @(r) r/2 + sigma_r^2/(2*dr^2);

A = zeros(M-1, M-1);
B = zeros(M-1, M-1);

for i=2:M-1
    A(i-1,i) = a(r_grid(i));
    A(i,i-1) = -b(r_grid(i));
    B(i-1,i) = -a(r_grid(i));
    B(i,i-1) = b(r_grid(i));
end

for i=1:M-1
    A(i,i) = 1/dt - c(r_grid(i));
    B(i,i) = 1/dt + c(r_grid(i));
end

% fill V by using % A*V(j) = B*V(j+1) + d;
[L,U] = lu(B);

for i=N+1:-1:2
    d = zeros(size(A,2),1); % what should d be?
    d(1) = b(r_grid(2))*(V(1,i-1) + V(1,i));
    d(end) = a(r_grid(end))*(V(end,i-1)+ V(end,i));
    V(2:M,i-1) = U\(L\(A*V(2:M,i) + d));
end

price = interp1(r_grid,V(:,1),r0);
end

%% function of bond price
function price = bond_price(T1, T2, sigma, k, theta, r)
    B= 1/k * (1 - exp(-k*(T2-T1)));
    A = exp((theta - sigma^2/(2*k^2))*(B - T2 + T1) - sigma^2/(4*k)*B^2);
    price = A * exp(-B*r);
end