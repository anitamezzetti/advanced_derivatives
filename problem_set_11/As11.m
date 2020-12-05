%% parameters
% for swaption
t = 0; % first date
r = 0.05; % term structure (flat)
Ts = 1; % first ex date
Tx = 2; % last ex date

% for underlying swap
tau = 0.25;
TN = 4; % last payment date
payment_dates = (0:16).*tau;

% corresponding forward rates
vols_fw = zeros(16,1);
vols_fw(2:7) = 0.2;
vols_fw(8:11) = .22;
vols_fw(12:16) = .24;

deltas = pi/2 * ((1:16) - 2)./14;