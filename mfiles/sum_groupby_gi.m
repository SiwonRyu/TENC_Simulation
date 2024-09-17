function [sum_g, index_gi_num] = sum_groupby_gi(index, vars)
Data    = [index,vars];
N_idx   = size(index,2);
K       = size(vars,2);
Data_sorted = sortrows(Data, 1:N_idx);

index_g = Data_sorted(:,1);
index_i = Data_sorted(:,2);
list_g = unique(index_g);
list_i = unique(index_i);

N_size  = numel(num2str(max(index_i)));
G_multiplier = power(10,N_size);
index_gi = G_multiplier*index_g+index_i; % gen. unique ID for (g,i)

[~, index_gi_num] = ismember(index_gi, unique(index_gi));
sum_g = [];
for k = 1:K
    sum_g = [sum_g,accumarray(index_gi_num,Data_sorted(:,N_idx+k))];
end


end