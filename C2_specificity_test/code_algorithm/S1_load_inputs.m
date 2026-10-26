clear; clc;

%% ===== Project paths =====
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);

lesion_rmap_dir = fullfile(project_dir, 'disease_data', 'Lesion_rmaps');
sensitivity_dir = fullfile(project_dir, 'disease_data', 'Lesion_sensitivitymap');
output_dir = fullfile(project_dir, 'output');
step_dir = fullfile(output_dir, 'S1_inputs');

if ~exist(step_dir, 'dir')
    mkdir(step_dir);
end

%% ===== Parameters =====
REQUIRE_SENSITIVITY_MAP = true;

%% ===== Identify the symptoms from the r-map subfolders =====
% One subfolder per symptom, and the folder name is the symptom name, so the
% number of symptoms follows the folder structure without any code change.
folder_listing = dir(lesion_rmap_dir);
symptom_names = {folder_listing([folder_listing.isdir]).name};
symptom_names = sort(symptom_names(~startsWith(symptom_names, '.')));
num_symptoms = numel(symptom_names);

if num_symptoms == 0
    error('No symptom subfolders found in: %s', lesion_rmap_dir);
end

fprintf('Found %d symptoms: %s\n', num_symptoms, strjoin(symptom_names, ', '));

%% ===== Load the lesion r-maps of every symptom =====
r_maps_per_symptom = cell(num_symptoms, 1);
names_per_symptom = cell(num_symptoms, 1);
symptom_per_lesion = cell(num_symptoms, 1);
image_geometry = [];

for symptom = 1:num_symptoms
    symptom_dir = fullfile(lesion_rmap_dir, symptom_names{symptom});
    [r_maps, lesion_names, shape, affine, header] = ...
        load_symptom_rmaps(symptom_dir);

    if isempty(image_geometry)
        image_geometry = struct('shape', shape, 'affine', affine, 'header', header);
    elseif ~isequal(shape, image_geometry.shape) || ...
            any(abs(affine(:) - image_geometry.affine(:)) > 1e-4)
        error(['Grid mismatch: %s has shape %s and an affine that differs ' ...
            'from the first symptom (%s).'], ...
            symptom_names{symptom}, mat2str(shape), mat2str(image_geometry.shape));
    end

    fprintf('  %-12s %3d lesion r-maps\n', symptom_names{symptom}, numel(lesion_names));
    r_maps_per_symptom{symptom} = r_maps;
    names_per_symptom{symptom} = lesion_names(:);
    symptom_per_lesion{symptom} = repmat(symptom, numel(lesion_names), 1);
end

lesion_r_maps = vertcat(r_maps_per_symptom{:});
lesion_names = vertcat(names_per_symptom{:});
lesion_symptom = vertcat(symptom_per_lesion{:});

num_lesions = size(lesion_r_maps, 1);
num_voxels = size(lesion_r_maps, 2);
fprintf('Total: %d lesions, %d voxels, grid %s\n', ...
    num_lesions, num_voxels, mat2str(image_geometry.shape));

%% ===== Load the rLNM sensitivity map of every symptom =====
sensitivity_maps = cell(1, num_symptoms);
sensitivity_files = cell(1, num_symptoms);
missing_symptoms = {};

for symptom = 1:num_symptoms
    map_file = find_sensitivity_map(sensitivity_dir, symptom_names{symptom});

    if isempty(map_file)
        missing_symptoms{end+1} = symptom_names{symptom}; %#ok<SAGROW>
        continue
    end

    sensitivity_files{symptom} = map_file;
    [map_data, map_affine] = io_nifti_read(map_file);

    if isequal(size(map_data), image_geometry.shape) && ...
            all(abs(map_affine(:) - image_geometry.affine(:)) < 1e-4)
        sampled_map = map_data;
        note = 'same grid';
    else
        sampled_map = resample_nearest(map_data, map_affine, ...
            image_geometry.affine, image_geometry.shape);
        note = 'resampled to the r-map grid';
    end

    sensitivity_vector = sampled_map(:)' > 0.5;
    sensitivity_maps{symptom} = sensitivity_vector;
    fprintf('  %-12s %d / %d voxels kept (%s)  <- %s\n', ...
        symptom_names{symptom}, sum(sensitivity_vector), num_voxels, note, map_file);
end

if ~isempty(missing_symptoms)
    report_missing_sensitivity(missing_symptoms, sensitivity_dir, REQUIRE_SENSITIVITY_MAP);
end

%% ===== Write the input summary and save =====
write_input_summary(fullfile(output_dir, 'S1_input_summary.txt'), ...
    lesion_rmap_dir, sensitivity_dir, symptom_names, lesion_symptom, lesion_names, ...
    sensitivity_maps, sensitivity_files, image_geometry);

save(fullfile(step_dir, 'S1_inputs.mat'), ...
    'lesion_r_maps', 'lesion_symptom', 'lesion_names', 'symptom_names', ...
    'sensitivity_maps', 'sensitivity_files', 'image_geometry', ...
    'num_symptoms', 'num_lesions', 'num_voxels');

fprintf('Saved: %s\n', fullfile(step_dir, 'S1_inputs.mat'));


%% ===== Local functions =====
function [r_maps, lesion_names, shape, affine, header] = load_symptom_rmaps(symptom_dir)
% Load every lesion r-map of one symptom into an (n_lesions x n_voxels) matrix.
% Accepts .nii and .nii.gz; when both exist for the same stem, .nii.gz wins.

file_listing = dir(fullfile(symptom_dir, '*.nii*'));
all_names = {file_listing(~[file_listing.isdir]).name};
keep = ~cellfun(@isempty, regexpi(all_names, '\.nii(\.gz)?$'));

if ~any(keep)
    error('No .nii/.nii.gz files in: %s', symptom_dir);
end
all_names = all_names(keep);

% A lesion can be stored as both x.nii and x.nii.gz. Keep only the .nii.gz
% copy, otherwise the same lesion is counted twice.
is_compressed = ~cellfun(@isempty, regexpi(all_names, '\.nii\.gz$'));
file_names = [all_names(is_compressed), all_names(~is_compressed)];
stems = regexprep(file_names, '\.nii(\.gz)?$', '', 'ignorecase');
[~, first_of_stem] = unique(stems, 'stable');
file_names = sort(file_names(first_of_stem));

if numel(all_names) > numel(file_names)
    warning(['%s: %d files found but only %d lesions after removing ' ...
        'duplicate .nii/.nii.gz pairs.'], symptom_dir, numel(all_names), numel(file_names));
end

[first_map, affine, header] = io_nifti_read(fullfile(symptom_dir, file_names{1}));
shape = size(first_map);
r_maps = zeros(numel(file_names), numel(first_map));
r_maps(1, :) = first_map(:)';

for lesion = 2:numel(file_names)
    [r_map, ~, ~] = io_nifti_read(fullfile(symptom_dir, file_names{lesion}));

    if ~isequal(size(r_map), shape)
        error('%s is %s but the first r-map was %s.', ...
            file_names{lesion}, mat2str(size(r_map)), mat2str(shape));
    end

    r_maps(lesion, :) = r_map(:)';
end

lesion_names = regexprep(file_names, '\.nii(\.gz)?$', '', 'ignorecase');
end


function map_file = find_sensitivity_map(sensitivity_dir, symptom_name)
% Find the rLNM sensitivity map of one symptom, case-insensitively.
%
% First a few conventional paths are tried, then all .nii.gz files below
% sensitivity_dir are scored by name (exact 100, prefix/suffix 80, contains 60).
% Only files whose name contains "significant" are candidates, so a continuous
% map such as true_overlap_8mm.nii.gz can never be picked by accident. A tie is
% an error rather than a silent pick.

map_file = '';
if ~exist(sensitivity_dir, 'dir')
    return
end

direct_candidates = { ...
    fullfile(sensitivity_dir, symptom_name, 'significant_8mm.nii.gz'), ...
    fullfile(sensitivity_dir, symptom_name, 'significant.nii.gz'), ...
    fullfile(sensitivity_dir, [symptom_name '.nii.gz'])};

for candidate = 1:numel(direct_candidates)
    if exist(direct_candidates{candidate}, 'file') == 2
        map_file = direct_candidates{candidate};
        return
    end
end

file_listing = dir(fullfile(sensitivity_dir, '**', '*.nii.gz'));
if isempty(file_listing)
    return
end

lower_name = lower(symptom_name);
scores = zeros(numel(file_listing), 1);

for file = 1:numel(file_listing)
    stem = regexprep(file_listing(file).name, '\.nii\.gz$', '', 'ignorecase');

    % Only significance masks are accepted here. Without this guard a missing
    % mask silently falls back to the continuous true_overlap map.
    if ~contains(lower(stem), 'significant')
        continue
    end

    [~, parent] = fileparts(file_listing(file).folder);
    scores(file) = max(name_match_score(stem, lower_name), ...
                       name_match_score(parent, lower_name));
end

best_score = max(scores);
if best_score <= 0
    return
end

best = find(scores == best_score);
if numel(best) > 1
    error(['Several sensitivity maps match "%s"; rename the wanted one to ' ...
        '<symptom>/significant_8mm.nii.gz.'], symptom_name);
end

map_file = fullfile(file_listing(best).folder, file_listing(best).name);
end


function score = name_match_score(token, lower_name)
% 0 = no match, 100 = identical, 80 = prefix or suffix, 60 = contains.
score = 0;
if isempty(token)
    return
end

lower_token = lower(token);
if strcmp(lower_token, lower_name)
    score = 100;
elseif startsWith(lower_token, lower_name) || endsWith(lower_token, lower_name)
    score = 80;
elseif contains(lower_token, lower_name)
    score = 60;
end
end


function report_missing_sensitivity(missing_symptoms, sensitivity_dir, require_map)
message = sprintf(['Missing sensitivity maps for: %s\n' ...
    '  Expected at: %s/<symptom>/significant_8mm.nii.gz'], ...
    strjoin(missing_symptoms, ', '), sensitivity_dir);

if require_map
    error('%s', message);
else
    warning('%s', message);
end
end


function write_input_summary(summary_file, lesion_rmap_dir, sensitivity_dir, ...
                             symptom_names, lesion_symptom, lesion_names, ...
                             sensitivity_maps, sensitivity_files, image_geometry)
file_id = fopen(summary_file, 'w');

fprintf(file_id, 'Selective specificity testing - input summary\n');
fprintf(file_id, '%s\n', repmat('=', 1, 72));
fprintf(file_id, 'Lesion r-maps   : %s\n', lesion_rmap_dir);
fprintf(file_id, 'Sensitivity maps: %s\n', sensitivity_dir);
fprintf(file_id, 'Grid            : %s\n\n', mat2str(image_geometry.shape));

fprintf(file_id, '%14s %9s  %s\n', 'Symptom', 'n_lesion', 'Sensitivity map');
for symptom = 1:numel(symptom_names)
    if isempty(sensitivity_maps{symptom})
        map_text = '(missing)';
    else
        map_text = sprintf('%d voxels  <-  %s', ...
            sum(sensitivity_maps{symptom}), sensitivity_files{symptom});
    end
    fprintf(file_id, '%14s %9d  %s\n', symptom_names{symptom}, ...
        sum(lesion_symptom == symptom), map_text);
end

fprintf(file_id, '\nTotal: %d lesions\n\n', numel(lesion_symptom));
fprintf(file_id, 'Lesion list (symptom, file):\n');
for lesion = 1:numel(lesion_names)
    fprintf(file_id, '  %-12s %s\n', symptom_names{lesion_symptom(lesion)}, ...
        lesion_names{lesion});
end

fclose(file_id);
end
