function [implied_vol,fval] = implied_volBS(S,K,r,T,q,C,sigma0)
%This function finds the implied volatility (slide 8 lecture 2)
%The implied vol is the vol I have to insert in the BS model in order to
%reproduce the observed price. Therefore it is the volatility for which the
%difference between the two prices is zero.
% C = observed prices
% implied_vol = implied vol
% fval = 

implied_vol = zeros(size(C));

%find the vol such that: call_price_BS - observed_call_price = 0

for i = 1:length(K)
    for j = 1:length(T)
        
        % function which must be zero: 
        funct = @(impl_vol) (CallBS(S,K(i),r,T(j),impl_vol,q) - C(i,j));
        % find the zero:
        [implied_vol(i,j),fval(i,j)] = fzero(funct,sigma0);
    end
end

if fval ~= zeros(size(fval))
    fprint("Error in finding the implied volatility")
end

end

