function [data, affine, header] = io_nifti_read(file_path)
% Read a NIfTI file.
%
% data    (x,y,z) double
% affine  4x4, 0-based voxel -> mm (matches nibabel's affine)
% header  pass back to io_nifti_write to reuse the affine
%
% Uses niftiread when available, and falls back to SPM.

if exist('niftiread', 'file') == 2 || exist('niftiread', 'builtin') == 5
    header = niftiinfo(file_path);
    data = double(squeeze(niftiread(header)));

    % Prefer the raw sform rows: they are 0-based, like nibabel's affine.
    if isfield(header, 'raw') && isfield(header.raw, 'srow_x')
        affine = double([header.raw.srow_x(:)'; ...
                         header.raw.srow_y(:)'; ...
                         header.raw.srow_z(:)'; 0 0 0 1]);
    else
        affine = double(header.Transform.T);
    end
    return
end

if exist('spm_vol', 'file') == 2 || exist('spm_vol', 'builtin') == 5
    volume = spm_vol(file_path);
    data = double(spm_read_vols(volume));
    affine = volume(1).mat;
    header = struct('spm', volume);
    return
end

error('Neither niftiread (Image Processing Toolbox) nor SPM is available.');
end
