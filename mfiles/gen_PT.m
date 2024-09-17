function [Data_individual, Data_dyadic] = gen_PT(N,G,p_D,theta,beta_true)
beta_true0 = beta_true(:,1);
beta_true1 = beta_true(:,2);

D = random('binomial', 1,p_D, N,G);
Sum_D = sum(D,1)-D;

% Generate Dyadic Array of Treatments
[i_idx, j_idx] = meshgrid(1:N,1:N);
Di_tmp = D(i_idx(:),:);
Dj_tmp = D(j_idx(:),:);
D_dyadic = [Di_tmp(:), Dj_tmp(:)];

% Generate Dyadic and Symmetric Error Term
eps1 = randn(N,N,G);
eps_tmp1 = eps1.*tril(ones(N),-1);
eps1 = eps_tmp1+permute(eps_tmp1,[2,1,3]);
eps1_adj = eps1 + 10000*eye(N); % To remove self-links

eps0 = randn(N,N,G);
eps_tmp0 = eps0.*tril(ones(N),-1);
eps0 = eps_tmp0+permute(eps_tmp0,[2,1,3]);
eps0_adj = eps0 + 10000*eye(N); % To remove self-links

[A1_pot, A1_obs, A0_obs] = gen_pot_link(N,G,theta,eps0,eps1,D);

% Generate Individual Error Term
u0 = randn(N,G);
u1 = randn(N,G);
Y1_obs          = gen_pot_outcome(beta_true1,A1_obs,eps1,u1,D);
Y0_obs          = gen_pot_outcome(beta_true0,A0_obs,eps0,u0,zeros(N,G));

% Pack    
Data_individual = combine_individual([D(:), Y0_obs(:), Y1_obs(:), Sum_D(:)],N,G);
Data_dyadic = remove_diagonal([D_dyadic, A0_obs(:), A1_obs(:)],N,G);
end




%% Subfunctions
function X_rev = remove_diagonal(X,N,G)
% input: X: (NxNxG) x K dyadic matrix
% output: X_rev: ((N-1)xNxG) x K

% Links in dyadic array
    [i_idx, j_idx, g_idx] = meshgrid(1:N,1:N,1:G);
    index_gij = [g_idx(:), i_idx(:), j_idx(:)];
    idx_dyadic = index_gij(:,2) ~= index_gij(:,3);
    X_tmp = [index_gij, X];
    X_rev = X_tmp(idx_dyadic,:);
end

function Y_observed = gen_pot_outcome(beta,A,eps,u,D)
    eps_i   = squeeze(sum(eps,2));
    q       = squeeze(sum(A.*permute(  D,[3,1,2]),2));
    r       = squeeze(sum(A.*permute(1-D,[3,1,2]),2));
    Y_observed = beta(1)+beta(2)*D+beta(3)*q+beta(4)*r + u+0.5*eps_i;
end

function Data_individual = combine_individual(X,N,G);
    [i_idx, j_idx] = meshgrid(1:G,1:N);
    Data_individual = [i_idx(:), j_idx(:), X];
end

function [A1_potential, A1_observed, A0_observed] = gen_pot_link(N,G,theta,eps0,eps1,D)
%----------- Input arguments ----------------------------------------------
% - d, e: 0 or 1 denoting counterfactual scenario of treatment assignments
% - eps: dyadic error term (NxNxG)
% - theta: coefficients for the binary response
%----------- Output -------------------------------------------------------
% - A_potential: binary array (NxNxG)

%%%%% This function generates potential links by a binary response model
% - A(d,e) = 1{theta(1) + theta(2)d + theta(3)e + theta(4)de >= eps}
% - P[A(d,e) = 1] = F(theta(1) + theta(2)d + theta(3)e + theta(4)de)

%%%%% Identical Links over all dyads
% If eps are identically distributed, then {A_ij(d,e)} are identical
% over all dyads. Specifically, it implies A_ij(d,e) ~ A_ji(d,e). However,
% this does not guarantee A_ij(d,e) = A_ji(d,e).

%%%%% Undirected Network
% The observed links are given by
% - A_ij =   A_ij(1,1)DiDj     + A_ij(0,0)(1-Di)(1-Dj)
%           +A_ij(1,0)Di(1-Dj) + A_ij(0,1)(1-Di)Dj
% - A_ji =   A_ji(1,1)DiDj     + A_ji(0,0)(1-Di)(1-Dj)
%           +A_ji(0,1)Di(1-Dj) + A_ji(1,0)(1-Di)Dj
% Therefore, for undirected link we need
% (a) A_ij(d,e) = A_ji(d,e) for d=e,
% (b) A_ij(1,0) = A_ji(0,1), and A_ij(0,1) = A_ji(1,0)
%
% (a) is satisfied when eps is symmetric: eps(i,j) = eps(j,i)
% (b) is satisfied when eps is symmetric and theta(2)=theta(3)

%%%%% Potential links might have an almost sure ordering for some theta.

    [i_idx, j_idx] = meshgrid(1:N,1:N);
    Di_tmp = D(i_idx(:),:);
    Dj_tmp = D(j_idx(:),:);
    
    K   = @(d,e) [1 d e d*e]*theta;
    h0  = @(Di,Dj) -1.5+0.3*Di+0.3*Dj-Di.*Dj;
    h1  = @(Di,Dj) h0(Di,Dj)-K(0,0);

    A0_observed = reshape( h0(Di_tmp, Dj_tmp) ,[N,N,G]) >= eps0;

    A1_tmp = @(d,e) K(d,e) + reshape( h1(Di_tmp, Dj_tmp) , [N,N,G] ) >= eps1;
    A1_potential = cat(4, A1_tmp(1,1), A1_tmp(1,0), A1_tmp(0,1), A1_tmp(0,0));
    
    DD = permute(D,[1,3,2]).*permute(D,[3,1,2]);
    ND = permute(1-D,[1,3,2]).*permute(D,[3,1,2]);
    DN = permute(D,[1,3,2]).*permute(1-D,[3,1,2]);
    NN = permute(1-D,[1,3,2]).*permute(1-D,[3,1,2]);

    A1_observed =   A1_potential(:,:,:,1).*DD + A1_potential(:,:,:,2).*DN ...
                +   A1_potential(:,:,:,2).*ND + A1_potential(:,:,:,4).*NN;

end