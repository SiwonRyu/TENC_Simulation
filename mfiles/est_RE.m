function [e_zeta,e_beta,e_pi,e_beta_X,e_beta_S] = est_RE(Data_dyad,Data_ind,N)
%%%%%%%%%%%% Estimation
% Sort
Data_dyad   = sortrows(Data_dyad, 1:3);
Data_ind    = sortrows(Data_ind, 1:2);

% Load from dyadic data
idx_d       = Data_dyad(:,1:3);
Di          = Data_dyad(:,4);
Dj          = Data_dyad(:,5);
A1          = Data_dyad(:,6);

% Load from individual data
idx_i       = Data_ind(:,1:2);
D           = Data_ind(:,3);
Y           = Data_ind(:,4);

% finite sample adjustment for standard error
L_ind       = size(idx_i,1); % # of observation of individual data GxN
L_dyad      = size(idx_d,1); % # of observation of dyadic data GxNx(N-1)
K           = 4; % # of regressors = 4 without covariates
qc_ind      = (L_ind-1)/(L_ind-K);
qc_dyad     = (L_dyad-1)/(L_dyad-K);

% Estimate zeta
W           = [ones(L_dyad,1), Di, Dj, Di.*Dj];
[e_zeta, ~, A1_hat, IF_zeta] = LS_dyadic_g_new(W,A1,idx_d(:,1));

% Compute q, r, q_obs, r_obs
[q_pred ,~] = sum_groupby_new(idx_d(:,1:2), A1_hat  .*   Dj );
[r_pred ,~] = sum_groupby_new(idx_d(:,1:2), A1_hat  .*(1-Dj));
[Dq     ,~] = sum_groupby_new(idx_d(:,1:2), W       .*   Dj );
[Dr     ,~] = sum_groupby_new(idx_d(:,1:2), W       .*(1-Dj));
[q_obs  ,~] = sum_groupby_new(idx_d(:,1:2), A1      .*   Dj );
[r_obs  ,~] = sum_groupby_new(idx_d(:,1:2), A1      .*(1-Dj));

Z           = [ones(L_ind,1), D, q_pred, r_pred];
Z_obs       = [ones(L_ind,1), D, q_obs , r_obs ];

% Compute, estimate the decomposition (pi)
[est_beta_1st,~,~,IF_beta]  = LS_dyadic_g_new(Z,Y,idx_i(:,1));
[e_beta,IF_beta_adj]        = se_adjust_RE(est_beta_1st(:,1),Z,IF_beta,Dq,Dr,IF_zeta);
[e_beta_X]                  = LS_dyadic_g_new(Z_obs,Y,idx_i(:,1));
[e_beta_S]                  = LS_dyadic_g_new(Z(:,1:2),Y,idx_i(:,1));

e_pi = compute_pi(qc_ind,N,e_beta(:,1),e_zeta(:,1),e_zeta(:,1),IF_beta_adj,IF_zeta,IF_zeta);
end