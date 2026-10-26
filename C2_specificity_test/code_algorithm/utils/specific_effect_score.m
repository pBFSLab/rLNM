function score = specific_effect_score(group_a, group_b)
% Voxel-wise specific effect: signed standardized mean difference between the
% target group and the contrast group, using a degrees-of-freedom weighted
% pooled standard deviation.
%
%   s_pooled = sqrt( ((n_a-1)*v_a + (n_b-1)*v_b) / (n_a+n_b-2) )
%   score    = (m_a - m_b) / s_pooled
%
% group_a  (n_a x V) target group
% group_b  (n_b x V) contrast group
% score    (1 x V), signed as target minus contrast

num_a = size(group_a, 1);
num_b = size(group_b, 1);

if num_a < 2 || num_b < 2
    error('Each group needs at least two lesions for the specific effect.');
end

mean_difference = mean(group_a, 1) - mean(group_b, 1);
pooled_variance = ((num_a - 1) * var(group_a, 0, 1) + ...
                   (num_b - 1) * var(group_b, 0, 1)) / (num_a + num_b - 2);
pooled_sd = sqrt(pooled_variance);

score = zeros(1, size(group_a, 2));
has_spread = pooled_sd > 0;
score(has_spread) = mean_difference(has_spread) ./ pooled_sd(has_spread);
end
