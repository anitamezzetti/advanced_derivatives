function [payoff] = black_scholes(S0,t,sigma,X,T,r,q)

d1 = (log(S0/X)+ (r - q + (sigma^2 / 2))* (T-t))/(sigma*sqrt(T-t));
d2 = d1 - sigma*sqrt(T-t);
payoff = S0 * exp(-q*(T-t))* normcdf(d1) - normcdf(d2) * exp(-r*(T-t)) * X;

end