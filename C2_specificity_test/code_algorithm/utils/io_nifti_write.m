function io_nifti_write(file_path, data, header)
% Write a NIfTI file in single precision, reusing the affine of the input image.
%
% file_path  .nii or .nii.gz (compressed automatically)
% data       (x,y,z) array
% header     third output of io_nifti_read (may be empty)

if endsWith(file_path, '.gz')
    file_path = file_path(1:end-3);   % niftiwrite adds the .gz itself
end
data = single(data);

if exist('niftiwrite', 'file') == 2 || exist('niftiwrite', 'builtin') == 5
    if isstruct(header) && isfield(header, 'raw')
        % The source image is often 4D (...,1) while we write 3D, so the
        % header dimensions have to be updated as well.
        output_header = header;
        data_size = size(data);

        output_header.ImageSize = data_size;
        output_header.PixelDimensions = output_header.PixelDimensions(1:numel(data_size));

        output_header.raw.dim(1) = numel(data_size);
        output_header.raw.dim(2:1+numel(data_size)) = data_size;
        output_header.raw.dim(2+numel(data_size):end) = 1;

        output_header.raw.pixdim(1) = 1;
        output_header.raw.pixdim(2:1+numel(data_size)) = output_header.PixelDimensions;
        output_header.raw.pixdim(2+numel(data_size):end) = 1;

        output_header.Datatype = 'single';
        niftiwrite(data, file_path, output_header, 'Compressed', true);
    else
        niftiwrite(data, file_path, 'Compressed', true);
    end
    return
end

if exist('spm_write_vol', 'file') == 2 || exist('spm_write_vol', 'builtin') == 5
    if ~isstruct(header) || ~isfield(header, 'spm')
        error('The SPM fallback needs the header returned by io_nifti_read.');
    end
    volume = header.spm(1);
    volume.fname = file_path;
    volume.dt = [16 0];
    spm_write_vol(volume, data);
    return
end

error('Neither niftiwrite nor SPM is available.');
end
