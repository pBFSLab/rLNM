function write_map(file_path, values, image_geometry, domain_index)
% Write a 1 x V vector as a NIfTI volume on the pipeline grid.
%
% values       1 x V (full grid) or 1 x numel(domain_index)
% domain_index optional linear indices into the full grid; omit or pass []
%              when values already covers the full grid.

num_voxels = prod(image_geometry.shape);

if nargin < 4 || isempty(domain_index)
    if numel(values) ~= num_voxels
        error('write_map: %d values for a %d-voxel grid.', numel(values), num_voxels);
    end
    full_vector = values(:)';
else
    if numel(values) ~= numel(domain_index)
        error('write_map: %d values for %d domain voxels.', ...
            numel(values), numel(domain_index));
    end
    full_vector = zeros(1, num_voxels);
    full_vector(domain_index) = values(:)';
end

io_nifti_write(file_path, reshape(full_vector, image_geometry.shape), ...
    image_geometry.header);
end
