function write_lesionwise_maps(root_dir, maps, lesion_names, symptom_names, ...
                               lesion_symptom, domain_index, image_geometry)
% Write one NIfTI file per lesion under <root_dir>/<symptom>/<lesion>.nii.gz.
% Files that already exist are kept.
num_written = 0;

for lesion = 1:size(maps, 1)
    symptom_dir = fullfile(root_dir, symptom_names{lesion_symptom(lesion)});
    file_path = fullfile(symptom_dir, [lesion_names{lesion} '.nii.gz']);

    if exist(file_path, 'file') == 2
        continue
    end
    if ~exist(symptom_dir, 'dir')
        mkdir(symptom_dir);
    end

    write_map(file_path, maps(lesion, :), image_geometry, domain_index);
    num_written = num_written + 1;
end

fprintf('  %d written, %d already present -> %s\n', ...
    num_written, size(maps, 1) - num_written, root_dir);
end
