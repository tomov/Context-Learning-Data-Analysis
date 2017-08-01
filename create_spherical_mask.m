function [sphere_mask, sphere_vol] = create_spherical_mask(x, y, z, r, filename)

% Create a mask from a list of ROI labels.
%
% INPUT:
% x, y, z = coordinates of the center voxel in MNI space
% r = radius in voxels i.e. the native coordinate space, NOT MNI space!
% filename = (optional) output filename where to save the .nii file. If not provided,
%            file won't be saved
%
% OUTPUT:
% sphere_mask = the resulting mask as a 3D binary vector
% sphere_vol = the corresponding SPM volume variable
%

% load whole-brain mask so we don't pick voxels outside of it
%
[mask, Vmask] = load_mask('masks/mask.nii');
Vmask.fname = ''; % !!!! IMPORTANT!!!!! b/c we return it
if exist('filename', 'var')
    Vmask.fname = filename;
end
sphere_vol = Vmask;

% convert sphere center to native coordinate space
%
cor = mni2cor([x y z], Vmask.mat);
x = cor(1);
y = cor(2);
z = cor(3);

% find boundary coordinates
%
[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);

% create the spherical mask
%
sphere_mask = zeros(size(mask));
for newx = floor(x - r) : ceil(x + r)
    if newx < min_x || newx > max_x, continue; end
    for newy = floor(y - r) : ceil(y + r)
        if newy < min_y || newy > max_y, continue; end
        for newz = floor(z - r) : ceil(z + r)
            if newz < min_z || newz > max_z, continue; end
            if ~mask(newx, newy, newz), continue; end
            if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > r^2, continue; end
            sphere_mask(newx, newy, newz) = 1;
        end
    end
end

sphere_mask = logical(sphere_mask);

% optionally save the mask
%
if exist('filename', 'var')
    sphere_vol.fname = filename;
    spm_write_vol(sphere_vol, sphere_mask);
end

