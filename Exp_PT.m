clc;clear all;close all;
format compact

seed = 2048;
pool = gcp(); % Get the current parallel pool
numWorkers = pool.NumWorkers;

addpath('mfiles');
% Set Parameters
N = 20;
p_D = 0.5; % Prob. of being treated

beta_true0 = [1;0;0.6;0.6];
beta_true1 = [2;1;0.8;0.6];
beta_true  = [beta_true0,beta_true1];  % Parameters in outcome response
Dbeta_true = beta_true1-[beta_true0(1);0;0;0];
theta_true = [-1;0.1;0.1;1];  % Parameters for link formation

% Comptue True Values
[mA1_true, mA0_true, HA_true, zeta0_true, zeta1_true, xi_true, pi_true] = compute_true_PT(Dbeta_true, theta_true, N);
%par_true= [zeta_true; beta_true; pi_true; beta_true];

G = 2000;
[Data_individual, Data_dyadic] = gen_PT(N,G,p_D,theta_true,beta_true);
[est_zeta1, est_xi, est_beta, est_beta_L, est_pi] = est_PT(Data_dyadic,Data_individual, N);

[zeta1_true, est_zeta1]
[xi_true, est_xi]
[Dbeta_true, est_beta]
[pi_true, est_pi]



%%




[zeta_true, est_zeta]
[xi_true, est_xi]

beta_true = [1;1;0.8;0.6];
[beta_true, est_beta]

