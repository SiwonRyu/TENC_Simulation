function [sumValues, idx_decom, uniqueIndex] = sum_groupby_new(index_new, V)

N_idx   = size(index_new,2);    % # of indices

% # of digits of the maximum value of index_i
N_digits = zeros(N_idx,1);

for i = 1:N_idx
    N_digits(i) = numel(num2str(max(index_new(:,i))));
end
N_digits_tmp = cumsum(flip(N_digits));
G_multiplier = flip(power(10, [0,N_digits_tmp']))';

index_comb = index_new*G_multiplier(2:end); % Gen. unique ID for (g,i)
index_set  = [index_comb index_new];


% 고유한 combined index 찾기
[uniqueIndex, ~, idx] = unique(index_comb);

% 각 combined index에 대해 세 번째 칼럼의 합계 구하기

sumValues = [];
for k = 1:size(V,2) % iteration over variables
    sumValues = [sumValues,accumarray(idx, V(:,k))];
end


idx_decom = [];
for i = 1:N_idx
    if i == 1
        idx_decom(:,i) = floor(uniqueIndex/G_multiplier(2));
    end
    if i > 1 & i < N_idx
        idx_decom(:,i) = floor(mod(uniqueIndex,G_multiplier(i))/G_multiplier(i+1));
    end       
    if i == N_idx
        idx_decom(:,i) = mod(uniqueIndex,G_multiplier(N_idx));
    end
end
end