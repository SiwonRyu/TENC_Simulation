function [est_zeta1, est_xi, est_beta, est_beta_L, est_pi] = est_PT(Data_dyadic,Data_individual,N)
%%%%%%%%%%%% Estimation Part

% Load from dyadic data
idx     = Data_dyadic(:,1:3);
Di      = Data_dyadic(:,4);
Dj      = Data_dyadic(:,5);
A0      = Data_dyadic(:,6);
A1      = Data_dyadic(:,7);
DA      = A1-A0;

% Load from individual data
idx_ind = Data_individual(:,[1,2]);
D       = Data_individual(:,3);
Y0      = Data_individual(:,4);
Y1      = Data_individual(:,5);
SD      = Data_individual(:,6);
DY      = Y1-Y0;

% Dyadic Regressions
W = [ones(size(Data_dyadic,1),1), Di, Dj, Di.*Dj];
[est_zeta1, ~, A1_hat, IF_zeta1]  = LS_dyadic_g(W,A1,idx);
[est_xi, ~, DA_hat, IF_xi]  = LS_dyadic_g(W,DA,idx);
A0_hat = A1_hat - DA_hat;
IF_zeta0 = IF_zeta1 - IF_xi;
%[est_zeta0, ~, A0_hat, IF_zeta0]  = LS_dyadic_g(W,A0,idx);
%est_xi = est_zeta1-est_zeta0;
%IF_xi = IF_zeta1-IF_zeta0;

% Compute H matrix
% H1 = [1 0 0 (N-1)*zeta(1);
%     0 1 0 (N-1)*zeta(2);
%     0 0 zeta(1)+zeta(3) -zeta(1);
%     0 0 zeta(2)+zeta(4) -zeta(2)];
xi0 = [0;0.020596;0.020596;0.17205];%
zeta0 = [0.066807;0.068859;0.068859;0.037439];%

A1_hat_true = W*zeta0 ;%
DA_hat_true = W*xi0;%
A0_hat_true = A1_hat_true - DA_hat_true;%

q1_hat   = sum_groupby_gi(Data_dyadic(:,1:2), A1_hat.*    Dj );
r1_hat   = sum_groupby_gi(Data_dyadic(:,1:2), A1_hat.*(1- Dj));
s0_hat   = sum_groupby_gi(Data_dyadic(:,1:2), A0_hat);

q1_hat_L = sum_groupby_gi(Data_dyadic(:,1:2), A1.*    Dj );
r1_hat_L = sum_groupby_gi(Data_dyadic(:,1:2), A1.*(1- Dj));
s0_hat_L = sum_groupby_gi(Data_dyadic(:,1:2), A0);

Z   = [ones(size(D,1),1), D, q1_hat  , r1_hat-s0_hat];
Z_L = [ones(size(D,1),1), D, q1_hat_L, r1_hat_L-s0_hat_L];

Diff_q_hat = sum_groupby_gi(idx, W.*Dj);
Diff_r_hat = sum_groupby_gi(idx, W.*(1-Dj));
Diff_s_hat = sum_groupby_gi(idx, W);

[est_beta, ~, ~, IF_beta]   = LS_ind_adj_g_PT(Z,DY,idx_ind,Diff_q_hat,Diff_r_hat,Diff_s_hat,IF_zeta1, IF_zeta0);
[est_beta_L, ~, ~, ~]       = LS_ind_g(Z_L,DY,idx_ind);

est_pi = Compute_pi(N, est_beta(:,1), est_zeta1(:,1), est_xi(:,1), IF_beta, IF_zeta1, IF_xi);

end

%% Subfunctions
function est_pi = Compute_pi(N, beta, zeta, xi, IF_beta, IF_zeta, IF_xi)
    G = size(IF_zeta,1);

    pi = [beta(2); 
          beta(4)*(N-1)*xi(2); 
          (beta(3)-beta(4))*zeta(1); 
          beta(3)*xi(3)];
    
    IF_DT = IF_beta(:,2);
    IF_DN = (IF_beta(:,4)*xi(2) + beta(4)*IF_xi(:,2))*(N-1);
    IF_IT = (IF_beta(:,3)-IF_beta(:,4))*zeta(1) + (beta(3)-beta(4))*IF_zeta(:,1);
    IF_IN = IF_beta(:,3)*xi(3) + beta(3)*IF_xi(:,3);
    IF = cat(2, IF_DT, IF_DN, IF_IT, IF_IN);

    Avar = IF'*IF/G;
    SE = sqrt(diag(Avar/G));
    t = pi./SE;
    p = 2*(1-cdf('T',abs(t),G-4)); % need to adjust DF

    crit  = abs(icdf('T',0.975,G-4));
    LB = pi - crit*SE;
    UB = pi + crit*SE;

    est_pi = [pi, SE, t, p, LB, UB];  
end