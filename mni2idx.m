function idx = mni2idx(mni,maskfile)
% Convert from MNI coordinates to index in the betas array. Passes through mni2cor first
%
% INPUT:
% mni = [x y z] voxel coordinates (also supports list of voxels, e.g. [x1 y1 z1; x2 y2 z2; ...])
% maskfile = .nii file of the betas mask, i.e. the one that was passed to ccnl_get_betas()
%
% OUTPUT:
% idx = index of voxel in betas array, such that voxel beta is betas(idx)
%

[mask, V] = load_mask(maskfile);

cor = mni2cor(mni, V.mat);

cor2idx = nan(size(mask));
cor2idx(mask) = 1:sum(mask(:));

i = sub2ind(size(cor2idx), cor(:,1), cor(:,2), cor(:,3));

idx = cor2idx(i);
