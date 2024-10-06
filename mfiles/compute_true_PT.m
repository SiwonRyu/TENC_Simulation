function [mA1, mA0, HA, zeta0, zeta1, xi, pi] = compute_true_PT(beta, theta, N)

lf_idx  = @(d,e) [1 d e d*e]*theta; % Index for linlf_idx formation
h0      = @(Di,Dj) -1.5+0.3*Di+0.3*Dj-Di.*Dj;
h1      = @(Di,Dj) h0(Di,Dj)-lf_idx(0,0);

mA0     = [normcdf(h0(1,1));
           normcdf(h0(1,0));
           normcdf(h0(0,1));
           normcdf(h0(0,0))];

mA1     = [normcdf(lf_idx(1,1)+h1(1,1));
           normcdf(lf_idx(1,0)+h1(1,0));
           normcdf(lf_idx(0,1)+h1(0,1));
           normcdf(lf_idx(0,0)+h1(0,0))];

HA      = mA1-mA0;

M       = [0 0 0 1; 0 1 0 -1; 0 0 1 -1; 1 -1 -1 1];
zeta0   = M*mA0;
zeta1   = M*mA1;
xi      = zeta1-zeta0;

% Type 1 decomposition
pi      = [beta(2); 
           beta(4)*(N-1)*xi(2); 
           (beta(3)-beta(4))*zeta1(1); 
           beta(3)*xi(3)];

% Type 2 decomposition
% pi = [ beta(2); 
%        beta(4)*(N-1)*xi(2); 
%        (beta(3)-beta(4))*(zeta(3)+zeta(1)); 
%        beta(4)*xi(3)];

end