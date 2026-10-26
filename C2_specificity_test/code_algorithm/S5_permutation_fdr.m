clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

output_dir = fullfile(project_dir, 'output');
input_file = fullfile(output_dir, 'S1_inputs', 'S1_inputs.mat');
domain_file = fullfile(output_dir, 'S2_symptom_maps', 'S2_symptom_maps.mat');
stats_file = fullfile(output_dir, 'S4_group_comparison', 'S4_group_comparison.mat');
residual_file = fullfile(output_dir, 'S3_lesion_specific', 'S3_lesion_specific.mat');
step_dir = fullfile(output_dir, 'S5_permutation_fdr');

if ~exist(stats_file, 'file')
    error('Run S4_group_comparison first. Missing: %s', stats_file);
end
if ~exist(step_dir, 'dir')
    mkdir(step_dir);
end

%% ===== Parameters =====
cfg = specificity_config();

%% ===== Load =====
load(stats_file, 't_stat', 'specific_effect', 'num_target');
load(residual_file, 'lesion_specific_maps');
load(input_file, 'symptom_names', 'lesion_symptom', 'image_geometry', ...
     'num_symptoms', 'num_voxels');
load(domain_file, 'domain_index');

num_lesions = numel(lesion_symptom);
num_domain_voxels = size(t_stat, 2);

%% ===== Permute the symptom labels =====
% Methods: "symptom labels were randomly reassigned across lesion-specific
% residual maps 10,000 times while preserving group sizes. For each permutation,
% the voxel-wise Welch's t-statistic ... was recalculated to generate a null
% distribution, against which the observed statistic was compared to obtain an
% empirical P value."
p_one = zeros(num_symptoms, num_domain_voxels);
p_two = zeros(num_symptoms, num_domain_voxels);

for symptom = 1:num_symptoms
    timer = tic;

    [p_one(symptom, :), p_two(symptom, :)] = label_permutation( ...
        lesion_specific_maps, num_target(symptom), t_stat(symptom, :), ...
        cfg.n_permutations, cfg.n_workers);

    fprintf('  %-12s %d permutations in %.0f s\n', ...
        symptom_names{symptom}, cfg.n_permutations, toc(timer));
end

%% ===== Benjamini-Hochberg FDR correction =====
% Methods: "Empirical P values were corrected across voxels using the
% Benjamini-Hochberg FDR procedure."
q_one = zeros(num_symptoms, num_domain_voxels);
q_two = zeros(num_symptoms, num_domain_voxels);

for symptom = 1:num_symptoms
    q_one(symptom, :) = bh_fdr(p_one(symptom, :));
    q_two(symptom, :) = bh_fdr(p_two(symptom, :));
end

if strcmpi(cfg.sided, 'one')
    p_values = p_one;
    q_values = q_one;
else
    p_values = p_two;
    q_values = q_two;
end

%% ===== Keep the significant voxels =====
% Methods: "Voxels with P_FDR < 0.05 were retained to define the
% symptom-specific map." The specific effect is reported by S4 but is not part
% of the retention rule.
significant = q_values < cfg.fdr_level;
specific_maps = zeros(size(specific_effect));
specific_maps(significant) = specific_effect(significant);

for symptom = 1:num_symptoms
    fprintf('  %-12s %d / %d voxels significant (q_min = %.4e)\n', ...
        symptom_names{symptom}, nnz(significant(symptom, :)), num_domain_voxels, ...
        min(q_values(symptom, :)));
end

%% ===== Write the maps =====
for symptom = 1:num_symptoms
    write_map(fullfile(output_dir, ...
        sprintf('S5_pvalue_%s.nii.gz', symptom_names{symptom})), ...
        p_values(symptom, :), image_geometry, domain_index);

    write_map(fullfile(output_dir, ...
        sprintf('S5_qvalue_%s.nii.gz', symptom_names{symptom})), ...
        q_values(symptom, :), image_geometry, domain_index);

    write_map(fullfile(output_dir, ...
        sprintf('S5_specific_map_%s.nii.gz', symptom_names{symptom})), ...
        specific_maps(symptom, :), image_geometry, domain_index);

    write_map(fullfile(output_dir, ...
        sprintf('S5_specific_mask_%s.nii.gz', symptom_names{symptom})), ...
        double(significant(symptom, :)), image_geometry, domain_index);
end

%% ===== Summary and save =====
write_summary(fullfile(output_dir, 'S5_summary.txt'), symptom_names, ...
    num_target, num_lesions, sum(significant, 2), min(q_values, [], 2), ...
    cfg.fdr_level, cfg.sided);

save(fullfile(step_dir, 'S5_permutation_fdr.mat'), ...
    'p_one', 'p_two', 'q_one', 'q_two', 'p_values', 'q_values', ...
    'significant', 'specific_maps');

fprintf('Saved: %s\n', fullfile(step_dir, 'S5_permutation_fdr.mat'));


%% ===== Local functions =====
function write_summary(summary_file, symptom_names, num_target, num_lesions, ...
                       num_significant, q_min, fdr_level, sided)
file_id = fopen(summary_file, 'w');

fprintf(file_id, 'Selective specificity testing - summary\n');
fprintf(file_id, '%s\n', repmat('=', 1, 72));
fprintf(file_id, 'Symptoms: %s\n', strjoin(symptom_names, ', '));
fprintf(file_id, 'Lesions : %d\n', num_lesions);
fprintf(file_id, 'FDR level: %g (%s-sided)\n', fdr_level, sided);
fprintf(file_id, 'Retention rule: q < %g\n', fdr_level);
fprintf(file_id, '\n%14s %9s %10s %12s\n', 'Symptom', 'n_target', 'n_sig', 'q_min');
for symptom = 1:numel(symptom_names)
    fprintf(file_id, '%14s %9d %10d %12.4e\n', symptom_names{symptom}, ...
        num_target(symptom), num_significant(symptom), q_min(symptom));
end
fprintf(file_id, '\nCore result: S5_specific_map_<symptom>.nii.gz\n');
fprintf(file_id, ['Changing the FDR level does not require re-running the ' ...
    'permutation; use specificity_rethreshold.\n']);

fclose(file_id);
end
