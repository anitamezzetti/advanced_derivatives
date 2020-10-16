function [C] = CallBS(S,K,r,T,sigma,q)
% CallBS finds the price of a Call option using BS
% We suppose t=0 => T-t = T
% See slides 19/20 Lecure 1

d1 = (log(S/K) + (r - q + sigma^2 / 2)*(T))/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);

C =  S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
end

