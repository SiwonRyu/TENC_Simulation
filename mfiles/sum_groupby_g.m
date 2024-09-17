function [sum_g, list_g] = sum_groupby_g(index, vars)
Data    = [index,vars];
N_idx   = size(index,2);
K       = size(vars,2);
Data_sorted = sortrows(Data, 1:N_idx);

index_g = Data_sorted(:,1);
list_g = unique(index_g);
[~, index_g_num] = ismember(index_g, list_g);

sum_g = [];
    for k = 1:K
        sum_g = [sum_g,accumarray(index_g_num,Data_sorted(:,N_idx+k))];
    end
end