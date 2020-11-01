%--------------------------------------------------------------------------
% Advanced Derivatives 
% Problem Set 5
% Marc-Aurèle Divernois
%
% Andreasen-Huge interpolation
%--------------------------------------------------------------------------
function outputC = AndreasenHuge(NoNaNindex,j,v,C,deltaT,deltaK,K)
    mids = round((NoNaNindex(1:end-1,1) + NoNaNindex(2:end,1)) / 2);
    LiGrid = length(K);
    theta = NaN(LiGrid,1);
    
    % Build piecewise function
    k = 1;
    for p = 1:length(mids)
        while k < mids(p)
            theta(k) = v(p);
            k = k+1;
        end
        theta(k:end) = v(end);
    end

    % Build matrix A
    z = 0.5 * deltaT/deltaK^2 * theta(2:end).^2;  %start at 2 because z_1 corresponds to second strike
    z(end) = [];                                  %remove z_n because not needed in matrix A
    transform = repmat(z,1,3).*[-1 2 -1];
    A = zeros(1,LiGrid);
    for m = 1:LiGrid-2
        A = [A; zeros(1,m-1) transform(m,:) zeros(1,LiGrid-2-m)];
    end
    A = [A; zeros(1,LiGrid)];
    A = A + eye(LiGrid);
    
    % Output next column --- method 1: invert A --- more precise and much
    % faster than method 2
    outputC = A\C(:,j);
    
    % Output next column --- method 2: solve A*C_t+1 = C_t    
    %fun2 = @(Cnext) (A*Cnext - C(:,j));
    %[outputC,fval2] = fsolve(fun2,C(:,j));
end