function [F, S] = ForwardRates(r, sigma, theta, K, dt, N_sim, N_x, N_n)
% Dynamics of Forward rates under Spot-LIBOR measure

% N_x: Number of exercise dates
% N_n: Number of swap payment dates

%dW_k = cos(theta_k) dW_1 + sin(theta_k) dW_2

F = zeros(N_sim, N_x+1, N_n);
S = zeros(N_sim, N_x+1);
F(:, 1, :) = r*ones(N_sim, 1, N_n);

for t=2:N_x+1
    temp = 0;
    for k=t:N_n
        dW_1 = randn(N_sim, 1);
        dW_2 = randn(N_sim, 1);
        
        dW_k = (cos(theta(k))*dW_1 + sin(theta(k))*dW_2)*sqrt(dt);
        
        % slide 19 Lecture 10
        rho_tk = cos(theta(t) - theta(k));
        
        % slide 6 Lecture 11
        temp = temp + dt*rho_tk*sigma(k)*F(:,t-1,k)./(1+dt*F(:,t-1,k));
        mu_k = sigma(k)*temp;
        F(:,t,k) = F(:,t-1,k).*exp((mu_k - 0.5*sigma(k)^2)*dt +...
           sigma(k)*dW_k);
    end
    %disp(F)
    
    % compute value of underlying swap (intrinsic val) at each ...
    % exercise date
    % Andersen algo: slide 9 Lecture 11
    P = zeros(N_sim, N_n);
    
    if (t > 1/dt)
        P(:, :) = cumprod(1./(1+dt.*F(:,t,:)),3);
        P(:,1:t-1) = 0;
        
        % Value of underlying swap at t (exercise date)
        S(:, t) = dt*diag(P * (reshape(K-F(:,t,:), [N_sim,N_n]))');
    end
end
end