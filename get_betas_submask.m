function betas = get_betas_submask(my_mask, betas)
% Given a custom mask and a bunch of beta vectors based on the subject group-level
% mask (i.e. get_betas('mask.nii', ...)), return beta vectors corresponding to the custom mask.
% This is in order to avoid having to use ccnl_get_beta for the custom
% mask and to do stuff like spotlight search.
%
% INPUT:
% my_mask = 3D binary mask whose betas we want to extract, in the space of
%           the subject group-level mask
% betas = [nTimePoints x nVoxels] matrix where nVoxels = # voxels in the
%         subject group-level mask (and they correspond to them too)
%
% OUTPUT:
% betas = [nTimePoints x mVoxels] matrix where mVoxels = # voxels in
%         my_mask

% Load the subject group-level mask
%
group_mask_filename = fullfile('masks', 'mask.nii'); % the subject group-level mask
[~, group_vol, group_mask] = load_mask(group_mask_filename);

% Make sure custom mask is in same coordinate space as the group-level mask
%
assert(isequal(size(my_mask), size(group_mask)));

% Get the indices of the voxels in the group-level mask; make sure they
% correspond to the voxels in the beta vector
%
group_inds = find(group_mask);
assert(numel(group_inds) == size(betas, 2));

% Get the indices of the voxels in the custom mask
%
my_inds = find(my_mask);

% For each voxel, check if it's part of the custom mask. If it is, then use
% its beta.
%
which = ismember(group_inds, my_inds);
betas = betas(:, which);

