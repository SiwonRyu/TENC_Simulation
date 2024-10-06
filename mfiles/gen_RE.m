function [Data_ind, Data_dyad] = gen_RE(N,G,p_D,theta,beta_true)
% Generate fake data for randomized experimental setting
D           = random('binomial', 1,p_D, N,G);
Degree      = sum(D,1)-D;

% Generate Dyadic Array of Treatments
[i, j]      = meshgrid(1:N,1:N);
Di_tmp      = D(i(:),:);
Dj_tmp      = D(j(:),:);
D_dyadic    = [Di_tmp(:), Dj_tmp(:)];

% Generate Dyadic and Symmetric Error Term
eps         = randn(N,N,G);
eps_tmp     = eps.*tril(ones(N),-1);
eps         = eps_tmp + permute(eps_tmp,[2,1,3]);
eps_adj     = eps + 10000*eye(N); % To remove self-links

% Generate Individual Error Term
u1          = randn(N,G);

% Potential Links (NxNxGx4)
[~,A1_obs]  = gen_pot_link(theta,eps_adj,D);
Y1_obs      = gen_outcome(beta_true,A1_obs,eps,u1,D);

% Pack    
Data_ind    = add_ind_index([D(:), Y1_obs(:), Degree(:)],N,G);
Data_dyad   = remove_diagonal([D_dyadic, A1_obs(:)],N,G);
end


%% Subfunctions
function X_rev = remove_diagonal(X,N,G)
% input     -X      : [(NxNxG) x K] matrix of dyadic variables
% output    -X_rev  : [((N-1)xNxG) x K] diagonal data removed
% Links in dyadic array
    [i,j,g] = meshgrid(1:N,1:N,1:G);
    index   = [g(:),i(:),j(:)];
    idx_del = index(:,2) ~= index(:,3);
    X_tmp   = [index, X];
    X_rev   = X_tmp(idx_del,:);
end
function Y_obs = gen_outcome(beta,A1,eps,u,D)
% Generate observed outcome according to the response function
    eps_i   = squeeze(sum(eps,2));
    q       = squeeze(sum(A1.*permute(  D,[3,1,2]),2));
    r       = squeeze(sum(A1.*permute(1-D,[3,1,2]),2));
    Y_obs   = beta(1) + beta(2)*D + beta(3)*q + beta(4)*r + u+0.5*eps_i;
end
function Data = add_ind_index(X,N,G)
% Attatch group-individual index (g,i) to individual-level data (X)
    [i,j]   = meshgrid(1:G,1:N);
    Data    = [i(:), j(:), X];
end

function [A_pot, A_obs] = gen_pot_link(theta,eps,D)
%----------- Input arguments ----------------------------------------------
% - d,e: 0 or 1 denoting counterfactual scenario of treatment assignments
% - eps: dyadic error term (NxNxG)
% - theta, omega: coefficients for the binary response
%----------- Output -------------------------------------------------------
% - A_pot: potential links of a binary array (NxNxG)
% - A_obs: observed links of a binary array (NxNxG)

%%%%% This function generates observed/potential links
% - I_b(d,e)    = [1, d, e, de]*b: a single index with coefficient b
% - A(d,e)      = 1{I_theta(d,e) >= eps}
% - P[A(d,e)=1] = F(I_theta(d,e))

%%%%% Identical Links over all dyads
% If eps are identically distributed, then {A_ij(d,e)} are identical over
% all dyads. It implies A_ij(d,e)~A_ji(d,e), but not A_ij(d,e)=A_ji(d,e).

%%%%% Undirected Network
% The observed links are given by
% - A_ij =   A_ij(1,1)DiDj     + A_ij(0,0)(1-Di)(1-Dj)
%           +A_ij(1,0)Di(1-Dj) + A_ij(0,1)(1-Di)Dj
% - A_ji =   A_ji(1,1)DiDj     + A_ji(0,0)(1-Di)(1-Dj)
%           +A_ji(0,1)Di(1-Dj) + A_ji(1,0)(1-Di)Dj
% Therefore, for undirected link, we need
% (a) A_ij(d,e) = A_ji(d,e) for d=e,
% (b) A_ij(1,0) = A_ji(0,1), and A_ij(0,1) = A_ji(1,0)
%
% (a) is satisfied when eps is symmetric: eps(i,j)=eps(j,i),
% (b) is satisfied when eps is symmetric and theta(2)=theta(3).

%%%%% Potential links might have an almost sure ordering for some theta
    
    I_b     = @(d,e,b) b(1) + b(2)*d + b(3)*e + b(4)*(d.*e);
    A_tmp   = @(d,e) I_b(d,e,theta) >=eps;

    % Potnetial link
    A_pot   = cat(4, A_tmp(1,1), A_tmp(1,0), A_tmp(0,1), A_tmp(0,0));
    
    DD      = permute(D,[1,3,2]).*permute(D,[3,1,2]);
    ND      = permute(1-D,[1,3,2]).*permute(D,[3,1,2]);
    DN      = permute(D,[1,3,2]).*permute(1-D,[3,1,2]);
    NN      = permute(1-D,[1,3,2]).*permute(1-D,[3,1,2]);

    % Bbserved link
    A_obs   = A_pot(:,:,:,1).*DD + A_pot(:,:,:,2).*DN ...
            + A_pot(:,:,:,2).*ND + A_pot(:,:,:,4).*NN;
end