clc;clear all;close all;

format compact

N = 20;
G = 100;
p_D = 0.3; % Prob. of being treated

theta = [-1;0.1;0.1;1];

D = random('binomial', 1, p_D, N,G);
%D = rand(N,G);


[i_idx, j_idx, g_idx] = meshgrid(1:N,1:N,1:G);
index_gij = [g_idx(:), i_idx(:), j_idx(:)];
idx_dyadic = index_gij(:,2) ~= index_gij(:,3);

[i_idx, j_idx] = meshgrid(1:N,1:N);
D_dyadic = [];
for g = 1:G
    D_dyadic = [D_dyadic; D(i_idx(:),g),D(j_idx(:),g)];
end
DD = D_dyadic(:,1).*D_dyadic(:,2);
DN = D_dyadic(:,1).*(1-D_dyadic(:,2));
ND = (1-D_dyadic(:,1)).*D_dyadic(:,2);
NN = (1-D_dyadic(:,1)).*(1-D_dyadic(:,2));

[eta0, eta_adj0] = gen_dyadic_error(N,G);
[eta1, eta_adj1] = gen_dyadic_error(N,G);

h0 = @(Di,Dj) -1.5+0.3*Di+0.3*Dj-Di.*Dj;
A0_obs = h0(D_dyadic(:,1), D_dyadic(:,2)) >= eta_adj0(:);

K   = @(d,e) [1 d e d*e]*theta;
h1  = @(Di,Dj) h0(Di,Dj)-K(0,0);

A1_pot_tmp = @(d,e) K(d,e) + h1(D_dyadic(:,1), D_dyadic(:,2)) >= eta_adj1(:);

A1_pot = cat(2, ...
    A1_pot_tmp(1,1), ...
    A1_pot_tmp(1,0), ...
    A1_pot_tmp(0,1), ...
    A1_pot_tmp(0,0));

A1_obs = A1_pot(:,1).*DD+A1_pot(:,2).*DN+A1_pot(:,3).*ND+A1_pot(:,4).*NN;

% Fake Dyadic Data
Data_dyadic = [index_gij, D_dyadic, A0_obs, A1_obs];
Data_dyadic = Data_dyadic(idx_dyadic,:);


% For step-by-step check
ATT_A1_0 = [normcdf(K(1,1)+h1(1,1))-normcdf(h0(1,1));
            normcdf(K(1,0)+h1(1,0))-normcdf(h0(1,0));
            normcdf(K(0,1)+h1(0,1))-normcdf(h0(0,1))];

ATT_A1_pot = [  sum((A1_pot(idx_dyadic,1)-A1_pot(idx_dyadic,4)).*DD(idx_dyadic))/sum(DD(idx_dyadic));
                sum((A1_pot(idx_dyadic,2)-A1_pot(idx_dyadic,4)).*DN(idx_dyadic))/sum(DN(idx_dyadic));
                sum((A1_pot(idx_dyadic,3)-A1_pot(idx_dyadic,4)).*ND(idx_dyadic))/sum(ND(idx_dyadic))];

EA1 = [ normcdf(K(1,1)+h1(1,1));
        normcdf(K(1,0)+h1(1,0));
        normcdf(K(0,1)+h1(0,1));
        normcdf(K(0,0)+h1(0,0))]


A1_obs_dyadic = Data_dyadic(:,7);
A0_obs_dyadic = Data_dyadic(:,6);
Dj_dyadic = Data_dyadic(:,5);






% Generate Individual Error Term
u0 = randn(N,G);
u1 = randn(N,G);

%eps1_   = squeeze(sum(eye(N).*eps1,2));
eta0_   = squeeze(sum(eta0,2));
eta1_   = squeeze(sum(eta1,2));

q1 = sum_group_by_gi(Data_dyadic(:,1:2), A1_obs_dyadic.*    Dj_dyadic );
r1 = sum_group_by_gi(Data_dyadic(:,1:2), A1_obs_dyadic.*(1- Dj_dyadic));
r0 = sum_group_by_gi(Data_dyadic(:,1:2), A0_obs_dyadic);

beta11 = 2;
beta10 = 1;
beta2 = 1;
beta3 = 0.8;
beta4 = 0.4;
beta = [beta11; beta10; beta2; beta3; beta4];

Y1_obs = beta11 + beta2*D(:)    + beta3*q1 + beta4*r1 + u1(:) + 0.5*eta1_(:);
Y0_obs = beta10                 + beta4*r0 + u0(:) + 0.5*eta0_(:);

index_gi = unique(Data_dyadic(:,1:2),'rows');
Data_ind = [index_gi, Y0_obs, Y1_obs, D(:)];




% Estimation
clearvars -except N G Data_dyadic ATT_A1_0 ATT_A1_pot EA1 Data_ind beta

% Load Dyadic Data of (GxIxJ)x{g,i,j,Di,Dj,Aij0,Aij1}
A1_obs_dyadic = Data_dyadic(:,7);
A0_obs_dyadic = Data_dyadic(:,6);
D_dyadic = Data_dyadic(:,4:5);
Dj_dyadic = Data_dyadic(:,5);



Di = Data_ind(:,5);
SDj = sum_group_by_gi(Data_dyadic(:,1:2), Dj_dyadic) - Di;
DiSDj = Di.*SDj;

Y1_obs = Data_ind(:,4);
Y0_obs = Data_ind(:,3);

DD = D_dyadic(:,1).*D_dyadic(:,2);
DN = D_dyadic(:,1).*(1-D_dyadic(:,2));
ND = (1-D_dyadic(:,1)).*D_dyadic(:,2);
NN = (1-D_dyadic(:,1)).*(1-D_dyadic(:,2));

% Tmp
m1_0 = mean(A1_obs_dyadic(logical(NN)));
SD1_DD = mean(A1_obs_dyadic(logical(DD)))-m1_0;
SD1_DN = mean(A1_obs_dyadic(logical(DN)))-m1_0;
SD1_ND = mean(A1_obs_dyadic(logical(ND)))-m1_0;

m0_0 = mean(A0_obs_dyadic(logical(NN)));
SD0_DD = mean(A0_obs_dyadic(logical(DD)))-m0_0;
SD0_DN = mean(A0_obs_dyadic(logical(DN)))-m0_0;
SD0_ND = mean(A0_obs_dyadic(logical(ND)))-m0_0;

ATT_A1_SD = [SD1_DD; SD1_DN; SD1_ND];% single diff
ATT_A1_DD = [SD1_DD - SD0_DD; SD1_DN - SD0_DN; SD1_ND - SD0_ND];% DiD


% Dyadic Regression
DA = A1_obs_dyadic - A0_obs_dyadic;
W = [ones(N*(N-1)*G,1) D_dyadic(:,1), D_dyadic(:,2), D_dyadic(:,1).*D_dyadic(:,2)];

xi = (W'*W)\(W'*DA);
zeta = (W'*W)\(W'*A1_obs_dyadic);
omega = (W'*W)\(W'*A0_obs_dyadic);
coefs = [xi zeta omega zeta-omega]

M1 = [0 1 1 1 ; 0 1 0 0 ; 0 0 1 0 ];
ATTs = [ATT_A1_0, ATT_A1_pot, ATT_A1_DD M1*xi] % compare: true/by pot/DD/xi



DA_hat = W*xi;
A1_hat = W*zeta;
A0_hat = W*omega;

H1 = [1 0 0 (N-1)*xi(1);
    0 1 0 (N-1)*xi(2);
    0 0 zeta(1)+zeta(3) xi(3)-zeta(1)-zeta(3);
    0 0 zeta(2)+zeta(4) xi(4)-zeta(2)-zeta(4)];




q1_hat = sum_group_by_gi(Data_dyadic(:,1:2), A1_hat.*    Dj_dyadic );
r1_hat = sum_group_by_gi(Data_dyadic(:,1:2), A1_hat.*(1- Dj_dyadic));
r0_hat = sum_group_by_gi(Data_dyadic(:,1:2), A0_hat);

q1_lim = sum_group_by_gi(Data_dyadic(:,1:2), A1_obs_dyadic.*    Dj_dyadic );
r1_lim = sum_group_by_gi(Data_dyadic(:,1:2), A1_obs_dyadic.*(1- Dj_dyadic));
r0_lim = sum_group_by_gi(Data_dyadic(:,1:2), A0_obs_dyadic);


DY = Y1_obs -Y0_obs;

Xtmp =[ones(N*G,1), Di, q1_hat, r1_hat-r0_hat];
beta_hat = (Xtmp'*Xtmp)\Xtmp'*DY;

Xlim =[ones(N*G,1), Di, q1_lim, r1_lim-r0_lim];
beta_lim = (Xlim'*Xlim)\Xlim'*DY;

ZW = [ones(N*G,1), Di, SDj, DiSDj];
delta_tmp = (ZW'*ZW)\(ZW'*DY)


beta_0 = [beta(1)-beta(2);beta(3);beta(4);beta(5)];
[beta_0, beta_hat, beta_lim, ]




function sum = sum_group_by_gi(index_dyadic, var_dyadic)
    
%[i_idx, g_idx] = meshgrid(1:N,1:G);
%index_gi = [g_idx(:), i_idx(:)]


    index_g = index_dyadic(:,1);
    index_i = index_dyadic(:,2);

    G = unique(index_g);
    N = unique(index_i);
    N_size  = numel(num2str(max(index_i)));
    G_multiplier = power(10,N_size);
    index_gi = G_multiplier*index_g+index_i; % gen. unique ID for (g,i)
    
    

    [~, index_gi_num] = ismember(index_gi, unique(index_gi));
    
    sum = accumarray(index_gi_num,var_dyadic);
end




function [eps, eps_adj] = gen_dyadic_error(N,G)
    eps = randn(N,N,G);
    %eps = (1/sqrt(2))*eps+permute((1/sqrt(2))*eps,[2,1,3]);
    eps_tmp = eps.*tril(ones(N),-1);
    eps = eps_tmp+permute(eps_tmp,[2,1,3]);
    eps_adj = eps + 10000*eye(N);
end
function [eps, eps_adj] = gen_dyadic_error_uni(N,G)
    eps = randn(N,N,G);
    eps_tmp = eps.*tril(ones(N),-1);
    eps = eps_tmp+permute(eps_tmp,[2,1,3]);
    eps_adj = eps + 10000*eye(N);
end

