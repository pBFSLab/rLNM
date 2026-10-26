function [p_one, p_two] = label_permutation(residual_maps, n_target, t_observed, ...
                                            n_perm, n_workers)
% Empirical P values from permuting symptom labels across lesions.
%
% Each permutation shuffles the lesion labels but keeps the two group sizes
% fixed, and recomputes the Welch t statistic:
%   p_two = (1 + #{|t_perm| >= |t_obs|}) / (1 + n_perm)   two-sided
%   p_one = (1 + #{ t_perm  >=  t_obs }) / (1 + n_perm)   one-sided
%
% residual_maps  (N x V) lesion-specific maps
% n_target       number of lesions in the target group
% t_observed     (1 x V) observed Welch t statistic
% n_perm         number of permutations; the loop runs exactly n_perm times
% n_workers      parallel workers; falls back to serial without the toolbox
%
% The permutations are split into n_workers chunks and chunk w starts its own
% random stream with rng(1000+w). Parallel and serial runs therefore draw the
% same permutation set, so the worker count only affects speed.

num_lesions = size(residual_maps, 1);
num_target = n_target;
num_contrast = num_lesions - n_target;
n_workers = max(1, round(n_workers));

% Split the permutations so the chunk counts add up to exactly n_perm.
chunk_counts = repmat(floor(n_perm / n_workers), 1, n_workers);
remainder = mod(n_perm, n_workers);
chunk_counts(1:remainder) = chunk_counts(1:remainder) + 1;
chunk_counts = chunk_counts(chunk_counts > 0);
num_chunks = numel(chunk_counts);

count_abs = zeros(num_chunks, size(residual_maps, 2));
count_ge = zeros(num_chunks, size(residual_maps, 2));

use_parallel = num_chunks > 1 && license('test', 'Distrib_Computing_Toolbox') ~= 0;
if use_parallel
    try
        if isempty(gcp('nocreate'))
            parpool('local', min(num_chunks, feature('numcores')), ...
                    'AutoAddClientPath', true);
        end
        parfor chunk = 1:num_chunks
            [count_abs(chunk, :), count_ge(chunk, :)] = permutation_chunk( ...
                residual_maps, num_target, num_contrast, t_observed, ...
                1000 + chunk, chunk_counts(chunk));
        end
    catch error_info
        warning('Parallel permutation failed; running serially.\n%s', ...
                getReport(error_info, 'basic'));
        use_parallel = false;
    end
end

if ~use_parallel
    for chunk = 1:num_chunks
        [count_abs(chunk, :), count_ge(chunk, :)] = permutation_chunk( ...
            residual_maps, num_target, num_contrast, t_observed, ...
            1000 + chunk, chunk_counts(chunk));
    end
end

p_two = (1 + sum(count_abs, 1)) / (1 + n_perm);
p_one = (1 + sum(count_ge, 1)) / (1 + n_perm);
end


function [count_abs, count_ge] = permutation_chunk(residual_maps, num_target, ...
                                                   num_contrast, t_observed, ...
                                                   seed, num_repeats)
% Run num_repeats permutations and accumulate the two counters per voxel.
num_voxels = size(residual_maps, 2);
count_abs = zeros(1, num_voxels);
count_ge = zeros(1, num_voxels);
abs_t_observed = abs(t_observed);

rng(seed, 'twister');
for repeat = 1:num_repeats
    order = randperm(num_target + num_contrast);
    group_a = residual_maps(order(1:num_target), :);
    group_b = residual_maps(order(num_target + 1:end), :);

    mean_difference = mean(group_a, 1) - mean(group_b, 1);
    squared_error = var(group_a, 0, 1) / num_target + ...
                    var(group_b, 0, 1) / num_contrast;

    t_perm = zeros(1, num_voxels);
    has_variance = squared_error > 0;
    t_perm(has_variance) = mean_difference(has_variance) ./ sqrt(squared_error(has_variance));

    count_abs = count_abs + (abs(t_perm) >= abs_t_observed);
    count_ge = count_ge + (t_perm >= t_observed);
end
end
