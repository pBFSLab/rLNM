clear; clc;

% Re-apply a different FDR level or specific-effect threshold to an existing
% result.
%
% The specific effect is stored by S4 and the permutation p/q values by S5, so
% changing a threshold never re-runs the permutation.

%% ===== Project paths =====
% This file lives in code_algorithm_latest/utils/, three levels down
% from the project root.
project_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

output_dir = fullfile(project_dir, 'output');
input_file = fullfile(output_dir, 'S1_inputs', 'S1_inputs.mat');
domain_file = fullfile(output_dir, 'S2_symptom_maps', 'S2_symptom_maps.mat');
stats_file = fullfile(output_dir, 'S4_group_comparison', 'S4_group_comparison.mat');
permutation_file = fullfile(output_dir, 'S5_permutation_fdr', 'S5_permutation_fdr.mat');
rethreshold_dir = fullfile(output_dir, 'S5_rethreshold');

%% ===== Parameters =====
cfg = specificity_config();
FDR_LEVEL = cfg.fdr_level;
SIDED = cfg.sided;

% Effect-size thresholds to scan. NaN means "no effect-size threshold".
SPECIFIC_EFFECT_THRESHOLDS = [0.5 0.8 1.0 1.2 NaN];

WRITE_MAPS = true;

%% ===== Load =====
load(stats_file, 'specific_effect');
load(permutation_file, 'q_one', 'q_two');
load(domain_file, 'domain_index');
load(input_file, 'symptom_names', 'image_geometry', 'num_symptoms');

if strcmpi(SIDED, 'one')
    q_values = q_one;
else
    q_values = q_two;
end

%% ===== Compare the thresholds =====
fprintf('FDR level %g, %s-sided\n\n', FDR_LEVEL, SIDED);
fprintf('%14s', 'Symptom');
for threshold = SPECIFIC_EFFECT_THRESHOLDS
    fprintf(' %10s', threshold_label(threshold));
end
fprintf('\n');

for symptom = 1:num_symptoms
    fprintf('%14s', symptom_names{symptom});
    for threshold = SPECIFIC_EFFECT_THRESHOLDS
        fprintf(' %10d', nnz(significant_voxels(q_values(symptom, :), ...
            specific_effect(symptom, :), FDR_LEVEL, threshold)));
    end
    fprintf('\n');
end

%% ===== Write the maps =====
if WRITE_MAPS
    if ~exist(rethreshold_dir, 'dir')
        mkdir(rethreshold_dir);
    end

    for threshold = SPECIFIC_EFFECT_THRESHOLDS
        tag = threshold_tag(threshold);

        for symptom = 1:num_symptoms
            is_significant = significant_voxels(q_values(symptom, :), ...
                specific_effect(symptom, :), FDR_LEVEL, threshold);

            specific_map = zeros(size(specific_effect(symptom, :)));
            specific_map(is_significant) = specific_effect(symptom, is_significant);

            write_map(fullfile(rethreshold_dir, ...
                sprintf('S5_specific_map%s_%s.nii.gz', tag, symptom_names{symptom})), ...
                specific_map, image_geometry, domain_index);

            write_map(fullfile(rethreshold_dir, ...
                sprintf('S5_specific_mask%s_%s.nii.gz', tag, symptom_names{symptom})), ...
                double(is_significant), image_geometry, domain_index);
        end
    end

    fprintf('\nWrote the thresholded maps to: %s\n', rethreshold_dir);
end


%% ===== Local functions =====
function is_significant = significant_voxels(q_values, specific_effect, fdr_level, threshold)
% NaN means "no effect-size threshold": keep every voxel that passes FDR.
if isnan(threshold)
    is_significant = q_values < fdr_level;
else
    is_significant = (q_values < fdr_level) & (specific_effect > threshold);
end
end


function tag = threshold_tag(threshold)
% File-name suffix, kept identical to S5_permutation_fdr.m.
if isnan(threshold)
    tag = '';
else
    tag = sprintf('_eff%g', threshold);
end
end


function label = threshold_label(threshold)
if isnan(threshold)
    label = 'FDR only';
else
    label = sprintf('effect>%g', threshold);
end
end
