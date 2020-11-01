%%% Advanced Derivatives 
%%% Black-Scholes formula
%%% Marc-Aurèle Divernois
function E = BlackScholes(V,t,sigma,D,T,r,q,type)
    d1 = (log(V/D)+ (r -q + sigma^2 / 2)* (T-t))/(sigma*sqrt(T-t));
    d2 = d1 - sigma*sqrt(T-t);
if type == "call"
    E = V * exp(-q*(T-t))* normcdf(d1) - normcdf(d2) * exp(-r*(T-t)) * D;
elseif type == "put"
    E = - V * exp(-q*(T-t))* normcdf(-d1) + normcdf(-d2) * exp(-r*(T-t)) * D;
else
    print("type error")
end
