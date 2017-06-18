function voxels = gen_rand_voxels(mask_filename, n_voxels)
% Generate random voxels within a given mask
%
% INPUT:
% mask_filename = path to .nii file with the mask
% n_voxels = how many voxels to generate
%
% OUTPUT:
% voxels = V x 3 list of [x y z] voxels in MNI space
%

[mask, Vmask, ~] = load_mask(mask_filename);

voxels = nan(n_voxels, 3);
for i = 1:n_voxels
    while true
        x = randi(size(mask, 1));
        y = randi(size(mask, 2));
        z = randi(size(mask, 3));
        if mask(x, y, z)
            voxels(i, :) = cor2mni([x y z], Vmask.mat);
            break;
        end
    end
end