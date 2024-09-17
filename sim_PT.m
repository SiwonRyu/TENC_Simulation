clc;clear all;close all;
format compact

seed = 2024;
pool = gcp(); % Get the current parallel pool
numWorkers = pool.NumWorkers;
filename = 'Result_PT0917';
sheetname = ['seed_' num2str(seed)];

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
par_true= [zeta1_true; xi_true; Dbeta_true; pi_true; Dbeta_true];

B       = 1000;
Glist   = [100 200 400 800 1600];
NG      = length(Glist);

par_est = zeros(20,6,NG,B);


tic(); gn = 0;
for G = Glist; gn = gn+1;
    parfor b = 1:B    
    % Create a unique random stream for each worker using an offset of the main seed
    stream = RandStream('mrg32k3a', 'Seed', seed + b);
    
    % Set the random stream for the current worker
    RandStream.setGlobalStream(stream);

    disp(['G =',num2str(G),', iteration: ',num2str(b),' / ',num2str(B), ' start'])    
    %zeta0 = [0.15866; 0.025405; 0.025405; 0.36979];
    
    [Data_individual, Data_dyadic] = gen_PT(N,G,p_D,theta_true,beta_true);
    [est_zeta, est_xi, est_beta, est_beta_L, est_pi] = est_PT(Data_dyadic,Data_individual, N);

    par_est(:,:,gn,b) = [est_zeta; est_xi; est_beta; est_pi; est_beta_L];
    end
end; 
clc; toc()

par_true'
[sum_par, sum_tot] = res_sum(par_est, par_true);
mean_res = sum_par(:,:,1)
med_res = sum_par(:,:,2);
mse_res = sum_par(:,:,3);
mae_res = sum_par(:,:,4);
bias_res = sum_par(:,:,5);
vari_res =sum_par(:,:,6);
cov_res =sum_par(:,:,7);

line_tmp = zeros(size(par_true,1)+1,1)';
mean_tmp = [Glist', mean_res];
mse_tmp = [Glist', mse_res];
cov_tmp = [Glist', cov_res];
res_tmp = [mean_tmp;
    0, par_true';
    line_tmp;
    mse_tmp; line_tmp;
    cov_tmp; line_tmp];
res_tmp = round(res_tmp,5);

T = array2table(res_tmp);
T.Properties.VariableNames = { 'G', ...
    'zeta1', 'zeta2', 'zeta3', 'zeta4', ...
    'xi1', 'xi2', 'xi3', 'xi4', ...
    'beta1', 'beta2', 'beta3', 'beta4', ...
    'pi1', 'pi2', 'pi3', 'pi4', ...
    'beta1L', 'beta2L', 'beta3L', 'beta4L'}

currentFolder = pwd;
writetable(T, [currentFolder '\Results\' filename '.xlsx'], 'Sheet', sheetname);



%%
%
% [sum_par, sum_tot] = res_sum(par_lim, mean(par_0,3));
% mean_res = sum_par(:,:,1)
% med_res = sum_par(:,:,2)
% mse_res = sum_par(:,:,3)
% mae_res = sum_par(:,:,4)
% bias_res = sum_par(:,:,5)
% vari_res =sum_par(:,:,6)
% toc()


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

% SS_bias_b  = sum(bias_b.^2,2);
% MAE_b_tot  = sum(MAE_b,2);
% MSE_b_tot  = permute(sum(mean(error_b.^2,4),1),[3,1,2]);
% var_b_tot  = sum(var_b,2);
% covr_b_min = min(covr_b,[],2);

% sum_par = cat(3,mean_b,med_b,MSE_b, MAE_b, bias_b, var_b, covr_b);
% sum_tot = [MSE_b_tot, MAE_b_tot, SS_bias_b, var_b_tot, covr_b_min];
sum_par = cat(3,mean_b,med_b,MSE_b, MAE_b, bias_b, var_b, covr_b);
sum_tot = 1;

end