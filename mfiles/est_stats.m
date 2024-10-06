function est_tab = est_stats(b,SE,level)
% Compute t-ratio, p-value, Confidence interval from coef and s.e.
% level: 1-significant level in [0,1]. Default = 95%
    
    if nargin == 2
        level = 0.95;
    end

    crit    = abs(icdf('Normal',1-(1-level)/2,0,1));
    t       = b./SE;
    %p      = 2*(1-cdf('T',abs(t),G-4)); % need to adjust DF
    %crit   = abs(icdf('T',0.975,G-4));
    p       = 2*(1-cdf('Normal',abs(t),0,1));
    
    LB      = b - crit*SE;
    UB      = b + crit*SE;
    
    est_tab = [b, SE, t, p, LB, UB];  
end