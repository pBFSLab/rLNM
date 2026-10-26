clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'output');
input_file = fullfile(output_dir, 'S1_inputs', 'S1_inputs.mat');
domain_file = fullfile(output_dir, 'S2_symptom_maps', 'S2_symptom_maps.mat');
residual_file = fullfile(output_dir, 'S3_lesion_specific', 'S3_lesion_specific.mat');
step_dir = fullfile(output_dir, 'S4_group_comparison');

if ~exist(residual_file, 'file')
    error('Run S3_lesion_specific first. Missing: %s', residual_file);
end
if ~exist(step_dir, 'dir')
    mkdir(step_dir);
end

%% ===== Load =====
load(residual_file, 'lesion_specific_maps');
load(input_file, 'lesion_symptom', 'symptom_names', 'image_geometry', ...
     'num_symptoms', 'num_voxels');
load(domain_file, 'domain_index');

num_lesions = numel(lesion_symptom);
num_domain_voxels = size(lesion_specific_maps, 2);

%% ===== Compare every symptom against all the others =====
% Methods: "For each SOI, the lesion-specific residual maps were divided into a
% target group ... and a contrast group ... The two groups were compared at each
% voxel using the Welch's t-test."
t_stat = zeros(num_symptoms, num_domain_voxels);
specific_effect = zeros(num_symptoms, num_domain_voxels);
num_target = zeros(num_symptoms, 1);

for symptom = 1:num_symptoms
    is_target = (lesion_symptom == symptom);
    num_target(symptom) = sum(is_target);

    fprintf('  %-12s target n=%d vs contrast n=%d\n', symptom_names{symptom}, ...
        num_target(symptom), num_lesions - num_target(symptom));

    t_stat(symptom, :) = welch_t(lesion_specific_maps(is_target, :), ...
                                 lesion_specific_maps(~is_target, :));
    specific_effect(symptom, :) = specific_effect_score( ...
        lesion_specific_maps(is_target, :), ...
        lesion_specific_maps(~is_target, :));
end

%% ===== Write the maps and save =====
for symptom = 1:num_symptoms
    write_map(fullfile(output_dir, ...
        sprintf('S4_tstat_%s.nii.gz', symptom_names{symptom})), ...
        t_stat(symptom, :), image_geometry, domain_index);

    write_map(fullfile(output_dir, ...
        sprintf('S4_specific_effect_%s.nii.gz', symptom_names{symptom})), ...
        specific_effect(symptom, :), image_geometry, domain_index);
end

save(fullfile(step_dir, 'S4_group_comparison.mat'), ...
    't_stat', 'specific_effect', 'num_target');

fprintf('Saved: %s\n', fullfile(step_dir, 'S4_group_comparison.mat'));
