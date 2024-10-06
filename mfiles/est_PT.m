function [e_zeta1,e_xi,e_beta,e_pi,e_beta_X,e_beta_S] = est_PT(Data_dyad,Data_ind,N)
%%%%%%%%%%%% Estimation Part

% Sort
Data_dyad   = sortrows(Data_dyad, 1:3);
Data_ind    = sortrows(Data_ind, 1:2);

% Load from dyadic data
idx_d       = Data_dyad(:,1:3);
Di          = Data_dyad(:,4);
Dj          = Data_dyad(:,5);
A0          = Data_dyad(:,6);
A1          = Data_dyad(:,7);
DA          = A1-A0;

% Load from individual data
idx_i       = Data_ind(:,[1,2]);
D           = Data_ind(:,3);
Y0          = Data_ind(:,4);
Y1          = Data_ind(:,5);
DY          = Y1-Y0;

% finite sample adjustment for standard error
L_ind       = size(idx_i,1); % # of observation of individual data GxN
L_dyad      = size(idx_d,1); % # of observation of dyadic data GxNx(N-1)
K           = 4; % # of regressors = 4 without covariates
qc_ind      = (L_ind-1)/(L_ind-K);
qc_dyad     = (L_dyad-1)/(L_dyad-K);

% Dyadic Regressions: estimate zeta, xi
W           = [ones(L_dyad,1), Di, Dj, Di.*Dj];
[e_zeta1  ,~, A1_hat, IF_zeta1]   = LS_dyadic_g_new(W,A1,idx_d(:,1));
[e_xi     ,~, DA_hat, IF_xi   ]   = LS_dyadic_g_new(W,DA,idx_d(:,1));
A0_hat      = A1_hat   - DA_hat;
IF_zeta0    = IF_zeta1 - IF_xi;

% Compute q, r, q_obs, r_obs
[q_pred ,~] = sum_groupby_new(idx_d(:,1:2), A1_hat  .*   Dj );
[r_pred ,~] = sum_groupby_new(idx_d(:,1:2), A1_hat  .*(1-Dj));
[s_pred ,~] = sum_groupby_new(idx_d(:,1:2), A0_hat          );
[Dq     ,~] = sum_groupby_new(idx_d(:,1:2), W       .*   Dj );
[Dr     ,~] = sum_groupby_new(idx_d(:,1:2), W       .*(1-Dj));
[Ds     ,~] = sum_groupby_new(idx_d(:,1:2), W               );
[q_obs  ,~] = sum_groupby_new(idx_d(:,1:2), A1      .*   Dj );
[r_obs  ,~] = sum_groupby_new(idx_d(:,1:2), A1      .*(1-Dj));
[s_obs  ,~] = sum_groupby_new(idx_d(:,1:2), A0              );

X           = [ones(size(D,1),1), D, q_pred, r_pred-s_pred];
X_obs       = [ones(size(D,1),1), D, q_obs , r_obs-s_obs  ];

% Compute, estimate the decomposition (pi)
[est_beta_1st,~, ~,IF_beta] = LS_dyadic_g_new(X,DY,idx_i(:,1));
[e_beta,IF_beta_adj]      = se_adjust_PT(est_beta_1st(:,1),X,IF_beta,Dq,Dr,Ds,IF_zeta0,IF_zeta1);
[e_beta_X,~,~,~]          = LS_dyadic_g_new(X_obs,DY,idx_i(:,1));
[e_beta_S,~,~,~]          = LS_dyadic_g_new(X(:,1:2),DY,idx_i(:,1));

e_pi = compute_pi(qc_ind,N,e_beta(:,1),e_zeta1(:,1),e_xi(:,1),IF_beta_adj,IF_zeta1,IF_xi);
end