function [payoff] = black_scholes(K,t,sigma,D,T,r,q)

d1 = (log(K/D)+ (r -q + sigma^2 / 2)* (T-t))/(sigma*sqrt(T-t));
d2 = d1 - sigma*sqrt(T-t);

payoff = K * exp(-q*(T-t))* normcdf(d1) - D* normcdf(d2) * exp(-r*(T-t));
end

