function [mA, HA, zeta, pi] = compute_true_RE(beta, theta, N)

lf_idx = @(d,e) [1 d e d*e]*theta; % Index for link formation
mA      = [normcdf(lf_idx(1,1));
           normcdf(lf_idx(1,0));
           normcdf(lf_idx(0,1));
           normcdf(lf_idx(0,0))];

HA      = [mA(4)
           mA(2)-mA(4)
           mA(3)-mA(4)
           mA(1)-mA(2)-mA(3)+mA(4)];

M       = [0 0 0 1; 0 1 0 -1; 0 0 1 -1; 1 -1 -1 1];

zeta    = M*mA;

% Type 1 decomposition
pi      = [beta(2); 
           beta(4)*(N-1)*zeta(2); 
           (beta(3)-beta(4))*zeta(1); 
           beta(3)*zeta(3)];

% Type 2 decomposition
% pi = [ beta(2); 
%        beta(4)*(N-1)*zeta(2); 
%        (beta(3)-beta(4))*(zeta(3)+zeta(1)); 
%        beta(4)*zeta(3)];
end