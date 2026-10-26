clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'output');
input_file = fullfile(output_dir, 'S1_inputs', 'S1_inputs.mat');
domain_file = fullfile(output_dir, 'S2_symptom_maps', 'S2_symptom_maps.mat');
step_dir = fullfile(output_dir, 'S3_lesion_specific');

if ~exist(domain_file, 'file')
    error('Run S2_symptom_maps first. Missing: %s', domain_file);
end
if ~exist(step_dir, 'dir')
    mkdir(step_dir);
end

%% ===== Load =====
load(input_file, 'lesion_r_maps', 'lesion_symptom', 'lesion_names', ...
     'symptom_names', 'image_geometry', 'num_symptoms');
load(domain_file, 'domain_index', 'num_domain_voxels');

lesion_maps = lesion_r_maps(:, domain_index);
num_lesions = size(lesion_maps, 1);

%% ===== Group sums for the leave-one-out averages =====
group_sum = zeros(num_symptoms, num_domain_voxels);
group_count = zeros(num_symptoms, 1);

for symptom = 1:num_symptoms
    in_group = (lesion_symptom == symptom);
    group_count(symptom) = sum(in_group);

    if group_count(symptom) < 2
        error(['Symptom %d has %d lesions; a leave-one-out average needs ' ...
            'at least two.'], symptom, group_count(symptom));
    end

    group_sum(symptom, :) = sum(lesion_maps(in_group, :), 1);
end

%% ===== LOLO regressor and OLS residual for every lesion =====
% Methods: "For each lesion, a leave-one-lesion-out (LOLO) shared-connectivity
% regressor was constructed ... The individual lesion connectivity map was then
% regressed on its corresponding LOLO shared-connectivity regressor across
% voxels using ordinary least squares."
lolo_regressor = zeros(num_lesions, num_domain_voxels);
lesion_specific_maps = zeros(num_lesions, num_domain_voxels);

for lesion = 1:num_lesions
    own_symptom = lesion_symptom(lesion);
    reference_maps = zeros(num_symptoms, num_domain_voxels);

    for symptom = 1:num_symptoms
        if symptom == own_symptom
            reference_maps(symptom, :) = (group_sum(symptom, :) - lesion_maps(lesion, :)) ...
                                         / (group_count(symptom) - 1);
        else
            reference_maps(symptom, :) = group_sum(symptom, :) / group_count(symptom);
        end
    end

    % x is left signed: the lesion map y is signed too, so the regression
    % decomposes y into shared and lesion-specific parts without clipping.
    x = mean(reference_maps, 1);
    lolo_regressor(lesion, :) = x;

    y = lesion_maps(lesion, :);
    x_centered = x - mean(x);
    y_centered = y - mean(y);
    sum_xx = sum(x_centered .^ 2);

    if sum_xx > 0
        slope = sum(x_centered .* y_centered) / sum_xx;
    else
        slope = 0;
    end

    lesion_specific_maps(lesion, :) = y_centered - slope * x_centered;
end

%% ===== Write the per-lesion maps and save =====
write_lesionwise_maps(fullfile(output_dir, 'S3_lolo'), lolo_regressor, ...
    lesion_names, symptom_names, lesion_symptom, domain_index, image_geometry);

write_lesionwise_maps(fullfile(output_dir, 'S3_lesion_specific'), ...
    lesion_specific_maps, lesion_names, symptom_names, lesion_symptom, ...
    domain_index, image_geometry);

save(fullfile(step_dir, 'S3_lesion_specific.mat'), ...
    'lolo_regressor', 'lesion_specific_maps');

fprintf('Saved: %s\n', fullfile(step_dir, 'S3_lesion_specific.mat'));
