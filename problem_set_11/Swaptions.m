function [CF] = Swaptions(S,H,D)
   % slide 10 to 12: Goal is to compute discounted cashflows as...
   %  D(0,T_e)*S(T_e, T_N)
   [~,num_dates] = size(S);
   CF_t = zeros(size(S));
   %disp(num_dates)
   for i=1:num_dates-1
       ex_idx = logical(S(:,end-i) > H);
       CF_t(ex_idx,end-i) = S(ex_idx,end-i).*D(ex_idx, end-i); 
       CF_t(ex_idx,end-i+1:end) = 0;
       for j=1:i
           ex_idx_next = logical(~(S(:,end-i+j-1) > H).*(S(:,end-i+j)>H));
           CF_t(ex_idx_next, end-i+j) = S(ex_idx_next,end-i+j).*...
               D(ex_idx_next,end-i+j);
       end
   end
    %disp(CF_t);
%     if(num_dates ==3)
%          disp(CF_t)
%     end
   CF = sum(CF_t,2);
   %disp(CF)
end