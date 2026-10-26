function q_values = bh_fdr(p_values)
% Benjamini-Hochberg step-up FDR correction.
%
%   q_(i) = min_{j >= i} ( m * p_(j) / j ),   with p sorted ascending
%
% Non-finite p values (NaN, Inf) come back as NaN.

p_shape = size(p_values);
p_values = double(p_values(:));
num_tests = numel(p_values);

is_bad = ~isfinite(p_values);
p_values(is_bad) = 1;                 % placeholder, restored to NaN below

[sorted_p, sort_order] = sort(p_values);

q_sorted = sorted_p .* (num_tests ./ (1:num_tests)');
q_sorted = flipud(cummin(flipud(q_sorted)));   % running minimum from the right
q_sorted = min(q_sorted, 1);

q_values = ones(num_tests, 1);
q_values(sort_order) = q_sorted;
q_values(is_bad) = NaN;
q_values = reshape(q_values, p_shape);
end
