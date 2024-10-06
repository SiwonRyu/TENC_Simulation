function est_pi = compute_pi(q, N, beta, zeta, xi, IF_beta, IF_zeta, IF_xi)
    pi = [ beta(2); 
           beta(4)*xi(2)*(N-1); 
           (beta(3)-beta(4))*zeta(1); 
           beta(3)*xi(3)];

    IF_DT = IF_beta(:,2);
    IF_DN = (IF_beta(:,4)*xi(2) + beta(4)*IF_xi(:,2))*(N-1);
    IF_IT = (IF_beta(:,3)-IF_beta(:,4))*zeta(1)+(beta(3)-beta(4))*IF_zeta(:,1);
    IF_IN = IF_beta(:,3)*xi(3) + beta(3)*IF_xi(:,3);
    IF    = cat(2, IF_DT, IF_DN, IF_IT, IF_IN);

    G = size(IF_zeta,1);
    Avar = q*IF'*IF/(G-1);
    SE = sqrt(diag(Avar/G));

    est_pi = est_stats(pi, SE);
end