function t_stat = welch_t(group_a, group_b)
% Voxel-wise Welch t statistic between two groups of lesion maps.
%
%   t = (m_a - m_b) / sqrt(v_a/n_a + v_b/n_b)     (equal variances not assumed)
%
% group_a  (n_a x V) target group
% group_b  (n_b x V) contrast group
% t_stat   (1 x V), signed as target minus contrast

num_a = size(group_a, 1);
num_b = size(group_b, 1);

if num_a < 2 || num_b < 2
    error('Each group needs at least two lesions for a Welch t-test.');
end

mean_difference = mean(group_a, 1) - mean(group_b, 1);
squared_error = var(group_a, 0, 1) / num_a + var(group_b, 0, 1) / num_b;

t_stat = zeros(1, size(group_a, 2));
has_variance = squared_error > 0;      % voxels with zero variance in both groups
t_stat(has_variance) = mean_difference(has_variance) ./ sqrt(squared_error(has_variance));
end
