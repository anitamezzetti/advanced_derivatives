%% Variables
S0 = 100;
T = 1;

t1 = 0.25; 
t2 = 0.5;
t3 = 0.75;
t4 = T;

ex_dates = [0, t1, t2, t3, t4];
dt = 0.25; % = delta t

K = 98;
r = 0;
q = 0.02; % dividend 

vol_c = 0.23; % constant volatility


%% stock (underlying) price matrix
N_paths = 1e5; % number of paths to simulate
S_mat = zeros(N_paths, 5); % N paths, 5 dates including t0 
S_mat(:, 1) = S0*ones(N_paths, 1); % starting prices are all 100
    
% prices with Brownian (slide 7 lecture 5):
for i=1:N_paths
    for j=2:5
        S_mat(i, j) = S_mat(i, j-1)*exp((r- q - vol_c^2/2)*dt +...
            vol_c*sqrt(dt)*normrnd(0,1));  
    end
end


%% Process A
A = zeros(N_paths, 4);
A(:,1) = S_mat(:, 2) - K;
A(:,2) = (S_mat(:, 2)+S_mat(:,3))/2 - K;
A(:,3) = (S_mat(:, 2)+S_mat(:,3) + S_mat(:,4))/3 - K;
A(:,4) = (S_mat(:, 2)+S_mat(:,3) + S_mat(:,4) + S_mat(:,5))/4 - K;


%% Cash flow matrix
C = zeros(N_paths, 4);
C(:,4) = max(A(:,4), zeros(N_paths, 1));

for i=3:-1:1
    intrinsic_vals = A(:, i); % exercise values
    ITM_prices_loc = find(intrinsic_vals > 0);
    
    ITM_prices = intrinsic_vals(ITM_prices_loc); % exercise values > 0
    discounted_cf = C(ITM_prices_loc, i+1)*exp(-r*dt);
    
    %X = [ones(length(ITM_prices)), ITM_prices, ITM_prices.^2,...
    %    ITM_prices.^3];
    
    % see https://ch.mathworks.com/help/matlab/data_analysis/linear-regression.html
    reg = polyfit(ITM_prices, discounted_cf, 3);
    fitted = polyval(reg, ITM_prices); % continuation values
    
    % Slide 11 lecture 5
    compare_matrix = zeros(N_paths,2);
    compare_matrix(ITM_prices_loc, 1) = ITM_prices;
    compare_matrix(ITM_prices_loc, 2) = fitted;
    
    % locations where exercise is more advantageous than continue
    ex_locs = find(compare_matrix(:, 1) > compare_matrix(:, 2));
    
    % Slide 12 and Slide 16 lecture 5
    C(ex_locs, i+1:end) = 0;
    C(ex_locs, i) = compare_matrix(ex_locs, 1);
    
end


price_t0 = mean(C(:,1)*exp(-r*dt));
disp(price_t0);


