function [par_est, par_true, time_sec] = simulation_RE(Glist, N, B, seed, p_D, beta, theta)
% Comptue True Values
[~,~,zeta, pi] = compute_true_RE(beta, theta, N);
par_true= [zeta; beta; pi; beta; 2; 1];

NG      = length(Glist);
par_est = zeros(18,6,NG,B);

tic(); gn = 0;
for G = Glist; gn = gn+1;
    parfor b = 1:B
    if ~isempty(seed)
        disp(['seed = ',num2str(seed), ' G =',num2str(G),', iteration: ',num2str(b),' / ',num2str(B), ' start'])    
        % Create a unique random stream for each worker
        stream = RandStream('mrg32k3a','Seed',seed + b);
        % Set random stream for the current worker
        RandStream.setGlobalStream(stream);
    else
        disp(['seed = none', ' G =',num2str(G),', iteration: ',num2str(b),' / ',num2str(B), ' start'])
    end

    [Data_ind, Data_dyad] = gen_RE(N,G,p_D,theta,beta)
    [est_zeta, est_beta, est_pi, est_beta_X, est_beta_S] ...
        = est_RE(Data_dyad,Data_ind,N);

    par_est(:,:,gn,b) = [est_zeta; est_beta; est_pi; est_beta_X; est_beta_S];
    end
end
clc; time_sec = toc()
end