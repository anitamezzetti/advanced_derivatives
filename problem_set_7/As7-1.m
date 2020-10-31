%% Variables
S0 = 100;
T = 1;

t1 = 0.25; 
t2 = 0.5;
t3 = 0.75;
t4 = T;

ex_dates = [t1, t2, t3, t4];
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
S_mat(:,2:end)=S_0 * cumprod(1 + (r-q)*dt + sqrt(dt)*vol_c.*randn(N_paths,4),2);

%% Process A
A = zeros(N_paths, 4);
A(:,1) = S_mat(:, 2);
A(:,2) = (S_mat(:, 2)+S_mat(:,3))/2;
A(:,3) = (S_mat(:, 2)+S_mat(:,3) + S_mat(:,4))/3;
A(:,4) = (S_mat(:, 2)+S_mat(:,3) + S_mat(:,4) + S_mat(:,5))/4;

I = max(A - K,0);
%% Cash flow matrix
C = zeros(N_paths, 4);
C(:,4) = max(A(:,4) - K, 0);

for i=3:-1:1
    intrinsic_vals = I(:,i); % exercise values
    ITM_prices_loc = find(intrinsic_vals > 0);
    
    ITM_prices = intrinsic_vals(ITM_prices_loc);
    
    S_to_regress = S_mat(ITM_prices_loc, i+1); 
    A_to_regress = A(ITM_prices_loc, i);
    discounted_cf = sum(exp(-r*(ex_dates(i+1:end) - ex_dates(i))).*...
        C(ITM_prices_loc, i+1:end),2);
    X = [ones(length(ITM_prices),1), S_to_regress, S_to_regress.^2,...
        S_to_regress.^3, A_to_regress, A_to_regress.^2, A_to_regress.^3];
    
    fitted = X * ((X'*X)\(X'*discounted_cf)); % continuation values
    
    % locations where exercise is more advantageous than continue
    ex_locs = find(ITM_prices > fitted);
    
    % Slide 12 and Slide 16 lecture 5
    C(ITM_prices_loc(ex_locs), i+1:end) = 0;
    C(ITM_prices_loc(ex_locs), i) = ITM_prices(ex_locs);
end


price_t0 = mean(sum(C.*exp(-r*ex_dates),2));
disp(price_t0);


