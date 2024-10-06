function [Data_ind, Data_dyad] = gen_PT(N,G,p_D,theta,omega,beta_true)
% Generate fake data for quasi-experimental setting under parallel trend
beta0       = beta_true(:,1);
beta1       = beta_true(:,2);

D           = random('binomial', 1,p_D, N,G);
Degree      = sum(D,1)-D;

% Generate Dyadic Array of Treatments
[i, j]      = meshgrid(1:N,1:N);
Di_tmp      = D(i(:),:);
Dj_tmp      = D(j(:),:);
D_dyadic    = [Di_tmp(:), Dj_tmp(:)];

% Generate Dyadic and Symmetric Error Term
eps1        = randn(N,N,G);
eps_tmp1    = eps1.*tril(ones(N),-1);
eps1        = eps_tmp1+permute(eps_tmp1,[2,1,3]);
eps1_adj    = eps1 + 10000*eye(N); % To remove self-links

eps0        = randn(N,N,G);
eps_tmp0    = eps0.*tril(ones(N),-1);
eps0        = eps_tmp0+permute(eps_tmp0,[2,1,3]);
eps0_adj    = eps0 + 10000*eye(N); % To remove self-links

[~, A1_obs, A0_obs] = gen_link(N,G,theta,omega,eps0_adj,eps1_adj,D);

% Generate Individual Error Term
u0          = randn(N,G);
u1          = randn(N,G);
Y1_obs      = gen_outcome(beta1,A1_obs,eps1,u1,D);
Y0_obs      = gen_outcome(beta0,A0_obs,eps0,u0,zeros(N,G));

% Pack    
Data_ind    = add_ind_index([D(:),Y0_obs(:),Y1_obs(:),Degree(:)],N,G);
Data_dyad   = remove_diagonal([D_dyadic,A0_obs(:),A1_obs(:)],N,G);
end


%% Subfunctions
function X_rev = remove_diagonal(X,N,G)
% input     -X      : [(NxNxG) x K] matrix of dyadic variables
% output:   -X_rev  : [((N-1)xNxG) x K] diagonal data removed
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

function [A1_pot, A1_obs, A0_obs] = gen_link(N,G,theta,omega,eps0,eps1,D)
%----------- Input arguments ----------------------------------------------
% - d,e: 0 or 1 denoting counterfactual scenario of treatment assignments
% - eps: dyadic error term (NxNxG)
% - theta, omega: coefficients for the binary response
%----------- Output -------------------------------------------------------
% - At_pot: potential links of a binary array (NxNxG)
% - At_obs: observed links of a binary array (NxNxG)

%%%%% This function generates observed/potential links
% - I_b(d,e)    = [1, d, e, de]*b: a single index with coefficient b
% - A0(d,e)     = A0    = 1{h0(Di,Dj) >= eps0}
% - A1(d,e)             = 1{h1(Di,Dj) + I_theta(d,e) >= eps0}
% - h0(d,e)     = I_omega(d,e)
% - h1(d,e)     = h0(d,e) - I_theta(0,0)
% - P[A0=1]     = F(I_omega(Di,Dj))
% - P[A1(d,e)=1]= F(I_omega(Di,Dj)-I_theta(0,0)+I_theta(d,e))

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
    
    % The distribution of (potential) links will be correlated with Di, Dj
    [i, j]  = meshgrid(1:N,1:N);
    Di_tmp  = D(i(:),:);
    Dj_tmp  = D(j(:),:);
    
    I_b     = @(d,e,b) b(1) + b(2)*d + b(3)*e + b(4)*(d.*e);
    h0      = @(d,e) I_b(d,e,omega);
    h1      = @(d,e) h0(d,e)-I_b(0,0,theta);

    % Indices for binary reseponse
    Sg_idx0 = reshape(h0(Di_tmp, Dj_tmp),[N,N,G]);
    Sg_idx1 = reshape(h1(Di_tmp, Dj_tmp),[N,N,G]);
    
    % Observed link at pre-treatment period
    A0_obs  = Sg_idx0 >= eps0;

    % Potnetial, and observed link at post-treatment period
    A1_tmp  = @(d,e) I_b(d,e,theta) + Sg_idx1 >= eps1;
    A1_pot  = cat(4,A1_tmp(1,1),A1_tmp(1,0),A1_tmp(0,1),A1_tmp(0,0));
    
    DD      = permute(D  ,[1,3,2]) .* permute(D  ,[3,1,2]);
    ND      = permute(1-D,[1,3,2]) .* permute(D  ,[3,1,2]);
    DN      = permute(D  ,[1,3,2]) .* permute(1-D,[3,1,2]);
    NN      = permute(1-D,[1,3,2]) .* permute(1-D,[3,1,2]);

    A1_obs  =   A1_pot(:,:,:,1).*DD + A1_pot(:,:,:,2).*DN ...
              + A1_pot(:,:,:,2).*ND + A1_pot(:,:,:,4).*NN;
end