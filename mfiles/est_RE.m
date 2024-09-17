function [est_zeta, est_beta, est_beta_L, est_pi] = est_RE(Data_dyadic,Data_individual,N)
%%%%%%%%%%%% Estimation Part

% Load from dyadic data
idx     = Data_dyadic(:,1:3);
Di      = Data_dyadic(:,4);
Dj      = Data_dyadic(:,5);
A1      = Data_dyadic(:,6);

% Load from individual data
idx_ind = Data_individual(:,[1,2]);
D       = Data_individual(:,3);
Y       = Data_individual(:,4);
SD      = Data_individual(:,5);

% Estimate zeta
% New version
W = [ones(size(Data_dyadic,1),1), Di, Dj, Di.*Dj];
[est_zeta, ~, A1_hat, IF_zeta] = LS_dyadic_g(W,A1,idx);

%zeta0 = [0.15866; 0.025405; 0.025405; 0.36979];
%A1_hat = W*zeta0;

% Compute H matrix
% H1 = [1 0 0 (N-1)*zeta(1);
%     0 1 0 (N-1)*zeta(2);
%     0 0 zeta(1)+zeta(3) -zeta(1);
%     0 0 zeta(2)+zeta(4) -zeta(2)];

% Compute q, r, q_obs, r_obs
q_hat   = sum_groupby_gi(Data_dyadic(:,1:2), A1_hat.*    Dj );
r_hat   = sum_groupby_gi(Data_dyadic(:,1:2), A1_hat.*(1- Dj));

Diff_q_hat = sum_groupby_gi(idx, W.*Dj);
Diff_r_hat = sum_groupby_gi(idx, W.*(1-Dj));

%q_hat_true   = sum_groupby_t(Data_dyadic(:,1:2), A1_hat_true.*    Dj ,1);
%r_hat_true   = sum_groupby_t(Data_dyadic(:,1:2), A1_hat_true.*(1- Dj),1);
q_hat_L = sum_groupby_gi(Data_dyadic(:,1:2), A1.*    Dj );
r_hat_L = sum_groupby_gi(Data_dyadic(:,1:2), A1.*(1- Dj));

Z   = [ones(size(D,1),1), D, q_hat, r_hat];
%Z   = [ones(N*G,1), D, q_hat_true, r_hat_true];
Z_L = [ones(size(D,1),1), D, q_hat_L, r_hat_L];

% [est_beta, ~, ~, IF_beta]   = LS(Z,Y);
% [est_beta_L, ~, ~, ~]       = LS(Z_L,Y);

% Compute, estimate the decomposition (pi)
% est_pi = Compute_pi(N, Data_dyadic, est_beta(:,1), est_zeta(:,1), IF_beta, IF_zeta);

%[est_beta, ~, ~, IF_beta]   = LS_ind_g(Z,Y,idx_ind);
[est_beta, ~, ~, IF_beta]   = LS_ind_adj_g(Z,Y,idx_ind,Diff_q_hat,Diff_r_hat,IF_zeta);
[est_beta_L, ~, ~, ~]       = LS_ind_g(Z_L,Y,idx_ind);

est_pi = Compute_pi(N, est_beta(:,1), est_zeta(:,1), IF_beta, IF_zeta);

%ZW = [ones(N*G,1), D, SD, D.*SD];
%delta_tmp = (ZW'*ZW)\(ZW'*Y)
%beta_est = [beta_true beta H1\delta_tmp delta_lim]

end

%% Subfunctions
function est_pi = Compute_pi(N, beta, zeta, IF_beta, IF_zeta)
    G = size(IF_zeta,1);

    pi = [beta(2); 
           beta(4)*(N-1)*zeta(2); 
           (beta(3)-beta(4))*zeta(1); 
           beta(3)*zeta(3)];

    IF_DT = IF_beta(:,2);
    IF_DN = (IF_beta(:,4)*zeta(2) + beta(4)*IF_zeta(:,2))*(N-1);
    IF_IT = (IF_beta(:,3)-IF_beta(:,4))*zeta(1) + (beta(3)-beta(4))*IF_zeta(:,1);
    IF_IN = IF_beta(:,3)*zeta(3) + beta(3)*IF_zeta(:,3);
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