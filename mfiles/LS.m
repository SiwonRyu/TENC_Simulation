function [Est_Table,r,Yhat,IF,Q] = LS(X,Y)
    L       = size(X,1); % Sum(N_g|g=1:G)
    K       = size(X,2); % Number of covariates
    b       = (X'*X)\X'*Y;
    Yhat    = X*b;
    r       = Y - X*b;
    Q       = (X'*X)/L;    
    Xr      = X.*r;
    q       = (L-1)/(L-K);
    IF      = (Q\Xr')';
    Avar    = q*IF'*IF/(L-1);
    SE      = sqrt(diag(Avar/L));
    Est_Table = est_stats(b, SE);
end