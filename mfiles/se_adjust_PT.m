function [Est_Table,IF_beta,R] = se_adjust_PT(b,X,IF,Dq,Dr,Ds,IF_zeta0,IF_zeta1)
%  Data = [index,X,Y];
%  Data_sorted = sortrows(Data, 1:2);
%  X_tmp = Data_sorted(:,3:7);
%  Y_tmp = Data_sorted

L       = size(X,1); % Number of obs (G x N x (N-1))
K       = size(X,2); % Number of covariates
q       = (L-1)/(L-K); % Finite Sample Adjustment

G       = size(IF,1);
R       = (X'*X)/G;

DZb1    = b(3).*Dq + b(4).*Dr;
DZb0    = b(4).*Ds;

% The ordering of Dq, Dr and X need to be the same    
IF_adj1  =  -(R\[( X'*DZb1 )/G]*IF_zeta1' )';
IF_adj0  =  -(R\[( X'*DZb0 )/G]*IF_zeta0' )';
IF_beta = IF + IF_adj1 - IF_adj0;

Avar    = q*IF_beta'*IF_beta/(G-1);
SE      = sqrt(diag(Avar/G));

Est_Table = est_stats(b, SE);
end







