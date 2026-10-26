% Run the whole specificity pipeline, S1 to S5, in order.
% Each step reads the .mat file written by the previous step, so the steps can
% also be run one by one. Helper functions live in utils/.

code_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(code_dir, 'utils'));

S1_load_inputs
S2_symptom_maps
S3_lesion_specific
S4_group_comparison
S5_permutation_fdr

fprintf('\nPipeline finished.\n');
