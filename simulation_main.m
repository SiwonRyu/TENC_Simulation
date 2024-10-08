clc;clear all;close all;
format compact

seed    = 2024;
filenm  = 'Result_1007_B10k_N10';
addpath('mfiles');

% Set paramteres for Design1
p_D     = 0.5; % Prob. of being treated
beta_R    = [2;1;0.8;0.6];  % Parameters in outcome response
theta_R   = [-1;0.1;0.1;1];  % Parameters for link formation

% Set paramteres for Design2
beta0_P = [1;0;0.6;0.6];
beta1_P = [2;1;0.8;0.6];
beta_P  = [beta0_P, beta1_P];  % Parameters in outcome response
theta_P = [-1;0.1;0.1;1];  % Parameters for link formation
omega_P = [-1.5; 0.3; 0.3;-1];

Glist   = [100 200 400 800 1600];
[par_est_R, par_true_R, time_sec_R] = simulation_RE(Glist,10,10000,seed,p_D,beta_R,theta_R);
[par_est_P, par_true_P, time_sec_P] = simulation_PT(Glist,10,10000,seed,p_D,beta_P,theta_P,omega_P);


% From Design 1
zeta_true_R   = par_true_R(1:4,:,:,:);
beta_true_R   = par_true_R(5:8,:,:,:);
pi_true_R     = par_true_R(9:12,:,:,:);
beta_S_true_R = par_true_R(17:18,:,:,:);

[sum_par_z_R, sum_tot_z_R]   = res_sum(par_est_R(1:4,:,:,:),   zeta_true_R);
[sum_par_b_R, sum_tot_b_R]   = res_sum(par_est_R(5:8,:,:,:),   beta_true_R);
[sum_par_p_R, sum_tot_p_R]   = res_sum(par_est_R(9:12,:,:,:),  pi_true_R);
[sum_par_bX_R, sum_tot_bX_R] = res_sum(par_est_R(13:16,:,:,:), beta_true_R);
[sum_par_bS_R, sum_tot_bS_R] = res_sum(par_est_R(17:18,:,:,:), beta_S_true_R);

mse_t_R = [sum_tot_z_R(:,1), sum_tot_b_R(:,1), sum_tot_p_R(:,1), sum_tot_bX_R(:,1), sum_tot_bS_R(:,1)];
mae_t_R = [sum_tot_z_R(:,2), sum_tot_b_R(:,2), sum_tot_p_R(:,2), sum_tot_bX_R(:,2), sum_tot_bS_R(:,2)];
ssb_t_R = [sum_tot_z_R(:,3), sum_tot_b_R(:,3), sum_tot_p_R(:,3), sum_tot_bX_R(:,3), sum_tot_bS_R(:,3)];
cov_t_R = [sum_tot_z_R(:,5), sum_tot_b_R(:,5), sum_tot_p_R(:,5), sum_tot_bX_R(:,5), sum_tot_bS_R(:,5)];

[sum_par_R, sum_tot_R] = res_sum(par_est_R, par_true_R);
mean_i_R  = sum_par_R(:,:,1);
med_i_R   = sum_par_R(:,:,2);
mse_i_R   = sum_par_R(:,:,3);
mae_i_R   = sum_par_R(:,:,4);
bias_i_R  = sum_par_R(:,:,5);
var_i_R   = sum_par_R(:,:,6);
cov_i_R   = sum_par_R(:,:,7);

% From Design 2
zeta1_true_P  = par_true_P(1:4,:,:,:);
xi_true_P     = par_true_P(5:8,:,:,:);
beta_true_P   = par_true_P(9:12,:,:,:);
pi_true_P     = par_true_P(13:16,:,:,:);
beta_S_true_P = par_true_P(21:22,:,:,:);

[sum_par_z_P,  sum_tot_z_P]  = res_sum(par_est_P(1:4,:,:,:), zeta1_true_P);
[sum_par_x_P,  sum_tot_x_P]  = res_sum(par_est_P(5:8,:,:,:), xi_true_P);
[sum_par_b_P,  sum_tot_b_P]  = res_sum(par_est_P(9:12,:,:,:), beta_true_P);
[sum_par_p_P,  sum_tot_p_P]  = res_sum(par_est_P(13:16,:,:,:), pi_true_P);
[sum_par_bL_P, sum_tot_bL_P] = res_sum(par_est_P(17:20,:,:,:), beta_true_P);
[sum_par_bS_P, sum_tot_bS_P] = res_sum(par_est_P(21:22,:,:,:), beta_S_true_P);

mse_t_P = [sum_tot_z_P(:,1), sum_tot_x_P(:,1), sum_tot_b_P(:,1), sum_tot_p_P(:,1), sum_tot_bL_P(:,1), sum_tot_bS_P(:,1)];
mae_t_P = [sum_tot_z_P(:,2), sum_tot_x_P(:,2), sum_tot_b_P(:,2), sum_tot_p_P(:,2), sum_tot_bL_P(:,2), sum_tot_bS_P(:,2)];
ssb_t_P = [sum_tot_z_P(:,3), sum_tot_x_P(:,3), sum_tot_b_P(:,3), sum_tot_p_P(:,3), sum_tot_bL_P(:,3), sum_tot_bS_P(:,3)];
cov_t_P = [sum_tot_z_P(:,5), sum_tot_x_P(:,5), sum_tot_b_P(:,5), sum_tot_p_P(:,5), sum_tot_bL_P(:,5), sum_tot_bS_P(:,5)];

[sum_par_P, sum_tot_P] = res_sum(par_est_P, par_true_P);
mean_i_P  = sum_par_P(:,:,1);
med_i_P   = sum_par_P(:,:,2);
mse_i_P   = sum_par_P(:,:,3);
mae_i_P   = sum_par_P(:,:,4);
bias_i_P  = sum_par_P(:,:,5);
var_i_P   = sum_par_P(:,:,6);
cov_i_P   = sum_par_P(:,:,7);

% Summarize
line    = zeros(size(par_true_R,1)+1+5,1)';
res_all_R = [ ...
    [Glist', mean_i_R, mae_t_R];
    0, par_true_R', zeros(1,5);
    line;
    [Glist', mse_i_R, mse_t_R];
    line;
    [Glist', cov_i_R, cov_t_R]];
res_all_R = round(res_all_R,5);

T_all_R = array2table(res_all_R);
T_all_R.Properties.VariableNames = { 'G', ...
    'zeta1', 'zeta2', 'zeta3', 'zeta4', ...
    'beta1', 'beta2', 'beta3', 'beta4', ...
    'pi1', 'pi2', 'pi3', 'pi4', ...
    'beta1L', 'beta2L', 'beta3L', 'beta4L', 'beta1S', 'beta2S', ...
    'zeta_T', 'beta_T', 'piT', 'bLT', 'bST'}

line = zeros(size(par_true_P,1)+1+6,1)';
res_all_P = [ ...
    [Glist', mean_i_P, mae_t_P];
    [0, par_true_P',zeros(1,6)];
    line;
    [Glist', mse_i_P, mse_t_P];
    line;
    [Glist', cov_i_P, cov_t_P]];
res_all_P = round(res_all_P,5);

T_all_P = array2table(res_all_P);
T_all_P.Properties.VariableNames = { 'G', ...
    'zeta1', 'zeta2', 'zeta3', 'zeta4', ...
    'xi1', 'xi2', 'xi3', 'xi4', ...
    'beta1', 'beta2', 'beta3', 'beta4', ...
    'pi1', 'pi2', 'pi3', 'pi4', ...
    'beta1L', 'beta2L', 'beta3L', 'beta4L', ...
    'beta1S', 'beta2S', ...
    'zeta_T', 'xi_T', 'beta_T', 'piT', 'bLT', 'bST'}

line = zeros(1,13);
res_zeta = [ ...
        [Glist', sum_par_z_R(:,:,1), sum_par_z_P(:,:,1), sum_par_x_P(:,:,1)];
        [0,      zeta_true_R',       zeta1_true_P',      xi_true_P'];
        line;
        [Glist', sum_par_z_R(:,:,2), sum_par_z_P(:,:,2), sum_par_x_P(:,:,2)];
        [0,      zeta_true_R',       zeta1_true_P',      xi_true_P'];
        line;
        [Glist', sum_par_z_R(:,:,3), sum_par_z_P(:,:,3), sum_par_x_P(:,:,3)];
        line;
        [Glist', sum_par_z_R(:,:,7), sum_par_z_P(:,:,7), sum_par_x_P(:,:,7)]];
res_zeta = round(res_zeta,5);
T_zeta = array2table(res_zeta);
T_zeta.Properties.VariableNames = { 'G', ...
    'zeta1R', 'zeta2R', 'zeta3R', 'zeta4R', ...
    'zeta1P', 'zeta2P', 'zeta3P', 'zeta4P', ...
    'xi1P'  , 'xi2P'  , 'xi3P'  , 'xi4P'}

line = zeros(1,9);
res_beta = [ ...
        [Glist', sum_par_b_R(:,:,1), sum_par_b_P(:,:,1)];
        [0,      beta_true_R',       beta_true_P'];
        line;
        [Glist', sum_par_b_R(:,:,2), sum_par_b_P(:,:,2)];
        [0,      beta_true_R',       beta_true_P'];
        line;
        line;
        [Glist', sum_par_b_R(:,:,3), sum_par_b_P(:,:,3)];
        line;
        [Glist', sum_par_b_R(:,:,7), sum_par_b_P(:,:,7)]];
res_beta = round(res_beta,5);
T_beta = array2table(res_beta);
T_beta.Properties.VariableNames = { 'G', ...
    'beta1R', 'beta2R', 'beta3R', 'beta4R', ...
    'beta1P', 'beta2P', 'beta3P', 'beta4P'}

line = zeros(1,9);
res_pi = [ ...
        [Glist', sum_par_p_R(:,:,1), sum_par_p_P(:,:,1)];
        [0,      pi_true_R',         pi_true_P'];
        line;
        [Glist', sum_par_p_R(:,:,2), sum_par_p_P(:,:,2)];
        [0,      pi_true_R',         pi_true_P'];
        line;
        [Glist', sum_par_p_R(:,:,3), sum_par_p_P(:,:,3)];
        line;
        [Glist', sum_par_p_R(:,:,7), sum_par_p_P(:,:,7)]];
res_pi = round(res_pi,5);
T_pi = array2table(res_pi);
T_pi.Properties.VariableNames = { 'G', ...
    'pi1_DTR', 'pi2_DNR', 'pi3_ITR', 'pi4_INR',...
    'pi1_DTP', 'pi2_DNP', 'pi3_ITP', 'pi4_INP'}

line = zeros(1,7);
res_compare = [ ...
        [Glist', mae_t_R(:,[2,4,5]), mae_t_P(:,[2,4,5])];
        line;
        [Glist', mse_t_R(:,[2,4,5]), mse_t_P(:,[2,4,5])];
        line;
        [Glist', cov_t_R(:,[2,4,5]), cov_t_P(:,[2,4,5])]]
res_compare = round(res_compare,5);
T_compare = array2table(res_compare);
T_compare.Properties.VariableNames = { 'G', ...
    'RR', 'XR', 'SR', ...
    'RP', 'XP', 'SP'}

cF = pwd;
writetable(T_all_R,     [cF '\Results\' filenm '.xlsx'], 'Sheet', 'AllR');
writetable(T_all_P,     [cF '\Results\' filenm '.xlsx'], 'Sheet', 'AllP');
writetable(T_pi,        [cF '\Results\' filenm '.xlsx'], 'Sheet', 'pi');
writetable(T_beta,      [cF '\Results\' filenm '.xlsx'], 'Sheet', 'beta');
writetable(T_zeta,      [cF '\Results\' filenm '.xlsx'], 'Sheet', 'zeta');
writetable(T_compare,   [cF '\Results\' filenm '.xlsx'], 'Sheet', 'compare');
%%
function [sum_par, sum_tot] = res_sum(res, par_true)
%%%%%%%%%%%% Summarize simulation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------- Input arguments ----------------------------------------------
%- res: container from simulation, M (# of parameters) x 4 (est/se/t/p) x G (# of G) x B(rep)
%- par_true: true parameter vector, M (# of parameters) x 1
%----------- Output -------------------------------------------------------
%- sum_par : mean/median/MSE/MAE/bias/variance/coverage
%- sum_tot : (agg.) MSE/MAE/Sum of squared bias/Sum of variance/coverage

Eb = @(res) permute(mean(res,4),[3,1,2]);
Mb = @(res) permute(median(res,4),[3,1,2]);

crit    = abs(icdf('Normal',0.975,0,1));
mean_b  = permute(mean(res(:,1,:,:),4),[3,1,2]);   % (Mx1xG) -> (GxM)
med_b   = permute(median(res(:,1,:,:),4),[3,1,2]); % (Mx1xG) -> (GxM)

mean_b  = Eb(res(:,1,:,:));
med_b   = Mb(res(:,1,:,:));

error_b = res(:,1,:,:) - par_true; % Mx1xGxB
bias_b  = permute(mean(error_b,4),[3,1,2]); 
MAE_b   = permute(mean(abs(error_b),4),[3,1,2]);
MSE_b   = permute(mean(error_b.^2,4),[3,1,2]);
dev_b   = res(:,1,:,:) - mean(res(:,1,:,:),4);
var_b   = permute(mean(dev_b.^2,4),[3,1,2]);

UB      = res(:,1,:,:) + crit*res(:,2,:,:);
LB      = res(:,1,:,:) - crit*res(:,2,:,:);
covr_b  = permute(mean((LB <= par_true).*(par_true <= UB), 4),[3,1,2]);

SS_bias_b  = sum(bias_b.^2,2);
MAE_b_tot  = mean(MAE_b,2);
MSE_b_tot  = permute(sum(mean(error_b.^2,4),1),[3,1,2]);
var_b_tot  = sum(var_b,2);
covr_b_min = permute(mean(prod((LB <= par_true).*(par_true <= UB),1), 4),[3,1,2]);

% sum_par = cat(3,mean_b,med_b,MSE_b, MAE_b, bias_b, var_b, covr_b);
sum_tot = [MSE_b_tot, MAE_b_tot, SS_bias_b, var_b_tot, covr_b_min];
sum_par = cat(3,mean_b,med_b,MSE_b, MAE_b, bias_b, var_b, covr_b);


end