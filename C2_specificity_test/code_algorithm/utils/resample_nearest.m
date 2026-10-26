function out = resample_nearest(src_data, src_affine, ref_affine, ref_shape, threshold)
% Resample src_data onto a reference grid with nearest-neighbour interpolation.
%
% Used to bring the 8 mm sensitivity maps onto the lesion r-map grid:
%   out(o) = src( inv(src_affine) * ref_affine * o ),  o = 0-based voxel index
%
% threshold  optional binarisation threshold (empty = keep the sampled values)

if nargin < 5
    threshold = [];
end

transform = double(src_affine) \ double(ref_affine);
src_size = size(src_data);
ref_shape = double(ref_shape(:))';

[ii, jj, kk] = ndgrid(0:ref_shape(1)-1, 0:ref_shape(2)-1, 0:ref_shape(3)-1);
x = transform(1,1)*ii + transform(1,2)*jj + transform(1,3)*kk + transform(1,4);
y = transform(2,1)*ii + transform(2,2)*jj + transform(2,3)*kk + transform(2,4);
z = transform(3,1)*ii + transform(3,2)*jj + transform(3,3)*kk + transform(3,4);

i = round(x) + 1;    % MATLAB indexing is 1-based
j = round(y) + 1;
k = round(z) + 1;

inside = (i >= 1) & (i <= src_size(1)) & ...
         (j >= 1) & (j <= src_size(2)) & ...
         (k >= 1) & (k <= src_size(3));

out = zeros(ref_shape);
linear_index = i(inside) + (j(inside) - 1) * src_size(1) + ...
               (k(inside) - 1) * src_size(1) * src_size(2);
out(inside) = src_data(linear_index);

if ~isempty(threshold)
    out = double(out > threshold);
end
end
