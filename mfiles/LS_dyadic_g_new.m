function [Est_Table,r,Yhat,IF,R] = LS_dyadic_g_new(X,Y,index)
L       = size(X,1); % Number of obs (G x N x (N-1))
K       = size(X,2); % Number of covariates
q       = (L-1)/(L-K); % Finite Sample Adjustment

b       = (X'*X)\X'*Y;
Yhat    = X*b;
r       = Y - X*b;

Xr      = sum_groupby_new(index, X.*r);
G           = size(Xr,1);
R         = (X'*X)/G;
IF        = (R\Xr')';

Avar_g    = q*IF'*IF/(G-1);
SE_g      = sqrt(diag(Avar_g/G));
Est_Table = est_stats(b, SE_g);
end