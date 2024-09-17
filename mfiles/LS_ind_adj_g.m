function [Est_Table,r,Yhat,IF,Q] = LS_ind_adj_g(X,Y,index,Dq,Dr,IF_zeta)
    Data = [index,X,Y];
    N_idx = size(index,2);
    Data_sorted = sortrows(Data, 1:N_idx);
    
    L       = size(X,1); % Sum(N_g|g=1:G)
    K       = size(X,2); % Number of covariates
    G       = max(unique(Data_sorted(:,1))); % Number of distinct groups

    % Mg:  G x 1 vector of # of dyads (without self link)
    N       = sum_groupby_g(index, ones(L,1));

    % Mg_long: (GxNx(N-1)) x 1 vector of # of dyads
    idx_g   = unique(index(:,1));
    N_long_tmp = [idx_g, N];
    N_long  = N_long_tmp(index(:,1),2);

    b       = (X'*X)\X'*Y;
    Yhat    = X*b;
    r       = Y - X*b;
    Q       = (X'*(X./N_long))/G;

    %Xr      = sum_groupby_t(index, X.*r, index_clustered);
    Xr      = sum_groupby_g(index, X.*r./N_long);
    IF      = (Q\Xr')';
    
    Diff_Zb = b(3).*Dq + b(4).*Dr;
    
    IF_adj  =  -(Q\[( X'*(Diff_Zb./N_long) )/G]*IF_zeta' )';
    IF_beta = IF + IF_adj;
    
    Avar    = IF_beta'*IF_beta/G;
    SE      = sqrt(diag(Avar/G));
    t       = b./SE;
    p       = 2*(1-cdf('T',abs(t),G-4)); % need to adjust DF
    
    crit        = abs(icdf('T',0.975,G-4));
    LB      = b - crit*SE;
    UB      = b + crit*SE;

    Est_Table = [b, SE, t, p, LB, UB];
end