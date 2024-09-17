clc;clear all;close all;
format compact

seed = 2048;
pool = gcp(); % Get the current parallel pool
numWorkers = pool.NumWorkers;

addpath('mfiles');

% Set Parameters
N = 20;
p_D = 0.5; % Prob. of being treated

beta_true = [2;1;0.8;0.6];  % Parameters in outcome response
theta_true = [-1;0.1;0.1;1];  % Parameters for link formation

% Comptue True Values
[mA_true, H_true, zeta_true, pi_true] = compute_true_RE(beta_true, theta_true, N);
par_true= [zeta_true; beta_true; pi_true; beta_true];

G = 200;

[Data_individual, Data_dyadic] = gen_RE(N,G,p_D,theta_true,beta_true);
[est_zeta1, est_beta1, est_beta_L1, est_pi1] = est_RE(Data_dyadic,Data_individual,N);

[zeta_true, est_zeta1]
[beta_true, est_beta1]
[pi_true, est_pi1]




