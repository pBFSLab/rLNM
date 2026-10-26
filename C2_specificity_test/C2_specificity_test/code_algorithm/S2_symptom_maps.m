clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'output');
input_file = fullfile(output_dir, 'S1_inputs', 'S1_inputs.mat');
step_dir = fullfile(output_dir, 'S2_symptom_maps');

if ~exist(input_file, 'file')
    error('Run S1_load_inputs first. Missing: %s', input_file);
end
if ~exist(step_dir, 'dir')
    mkdir(step_dir);
end

%% ===== Parameters =====
% The analysis domain is the brain mask only, so the voxels entering the FDR
% correction are fixed and are not selected by their own data.
BRAIN_MASK_FILE = fullfile(project_dir, 'template', 'norm_8mm.nii.gz');
if ~exist(BRAIN_MASK_FILE, 'file')
    % Sibling rLNM sensitivity-analysis template as a fallback.
    BRAIN_MASK_FILE = fullfile(project_dir, '..', 'S1_sensitivity_test', ...
                               'template', 'norm_8mm.nii.gz');
end

%% ===== Load =====
load(input_file, 'lesion_r_maps', 'lesion_symptom', 'symptom_names', ...
     'sensitivity_maps', 'image_geometry', 'num_symptoms', 'num_voxels');

%% ===== Symptom-level maps =====
% Methods: "For each SOI, the corresponding lesion-wise maps were averaged to
% derive a voxel-wise symptom-level connectivity map."
symptom_level_maps = zeros(num_symptoms, num_voxels);

for symptom = 1:num_symptoms
    symptom_level_maps(symptom, :) = ...
        mean(lesion_r_maps(lesion_symptom == symptom, :), 1);
end

%% ===== Positive restriction =====
% Methods: "The resulting symptom-level maps were then restricted to positive
% connectivity ..."
symptom_positive_maps = symptom_level_maps;
symptom_positive_maps(symptom_positive_maps <= 0) = 0;

for symptom = 1:num_symptoms
    fprintf('  %-12s %d / %d voxels positive\n', symptom_names{symptom}, ...
        nnz(symptom_positive_maps(symptom, :)), num_voxels);
end

%% ===== Analysis domain = brain mask =====
if isempty(BRAIN_MASK_FILE) || exist(BRAIN_MASK_FILE, 'file') ~= 2
    error('Brain mask not found: %s', BRAIN_MASK_FILE);
end

brain_mask = io_nifti_read(BRAIN_MASK_FILE);
brain_mask = logical(reshape(brain_mask, 1, []));

if numel(brain_mask) ~= num_voxels
    error('The brain mask has %d voxels but the r-maps have %d.', ...
        numel(brain_mask), num_voxels);
end

analysis_domain = brain_mask;
domain_index = find(analysis_domain);
num_domain_voxels = numel(domain_index);

if num_domain_voxels == 0
    error('The analysis domain is empty.');
end

fprintf('Analysis domain (brain mask): %d / %d voxels <- %s\n', ...
    num_domain_voxels, num_voxels, BRAIN_MASK_FILE);

%% ===== Cross-symptom map =====
% Methods: "... and averaged voxel-wise with equal weight across symptoms.
% The averaged map was then masked by the intersection of the binary
% significance masks obtained from null-model sensitivity testing for all
% symptoms in the comparison."
cross_symptom_map = mean(symptom_positive_maps, 1);

sensitivity_intersection = true(1, num_voxels);
for symptom = 1:num_symptoms
    if ~isempty(sensitivity_maps{symptom})
        sensitivity_intersection = ...
            sensitivity_intersection & sensitivity_maps{symptom};
    end
end

cross_symptom_map(~sensitivity_intersection) = 0;

fprintf('Sensitivity intersection: %d voxels; cross-symptom map: %d non-zero voxels\n', ...
    nnz(sensitivity_intersection), nnz(cross_symptom_map));

%% ===== Write the maps and save =====
for symptom = 1:num_symptoms
    write_map(fullfile(output_dir, ...
        sprintf('S2_symptom_map_%s.nii.gz', symptom_names{symptom})), ...
        symptom_level_maps(symptom, :), image_geometry);

    write_map(fullfile(output_dir, ...
        sprintf('S2_symptom_map_positive_%s.nii.gz', symptom_names{symptom})), ...
        symptom_positive_maps(symptom, :), image_geometry);
end

write_map(fullfile(output_dir, 'S2_domain_mask.nii.gz'), ...
    double(analysis_domain), image_geometry);
write_map(fullfile(output_dir, 'S2_cross_symptom_map.nii.gz'), ...
    cross_symptom_map, image_geometry);
write_map(fullfile(output_dir, 'S2_sens_intersection.nii.gz'), ...
    double(sensitivity_intersection), image_geometry);

save(fullfile(step_dir, 'S2_symptom_maps.mat'), ...
    'symptom_level_maps', 'symptom_positive_maps', 'analysis_domain', ...
    'domain_index', 'num_domain_voxels', 'sensitivity_intersection', ...
    'cross_symptom_map');

fprintf('Saved: %s\n', fullfile(step_dir, 'S2_symptom_maps.mat'));
