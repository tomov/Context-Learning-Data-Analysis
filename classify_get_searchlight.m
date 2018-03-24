function Searchlight = classify_get_searchlight(data, metadata, which_rows, event, x, y, z, r, use_pregenerated_activations, use_tmaps, use_nosmooth)

% Get searclight betas for classifier. Similar to classify_get_searchlights.m except for single searchlight.
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
% event = feedback_onset or trial_onset
% x, y, z = vectors of voxel coordinates in group-level mask coordinate space
% r = sphere radius in group-level coordinate space (i.e. voxels, NOT mm!)
% use_pregenerated_activations = whether to use the pregenerated whole-brain
%     betas (or t-values) to extract the submasks (faster but requires lots of memory to
%     preload them and fails on NCF...)
%     setting to false is SUPER SLOW...
% use_tmaps = if true, use tmaps; otherwise, use betas for neural data
% use_nosmooth = whether to use the non-smoothed neural data
%
% OUTPUT:
% Searchlight = searchlight betas 

searchlight_idx = 0;

assert(ismember(event, {'trial_onset', 'feedback_onset'}));

disp('Computing searchlights...');
tic

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end

[mask, Vmask] = load_mask('masks/mask.nii');
Vmask.fname = 'masks/searchlight.nii'; % !!!! IMPORTANT! in case we save it

[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);


if use_pregenerated_activations
    % Load neural data for given event. Warning -- takes up ~8.5 G and
    % NCF doens't like that
    %
    whole_brain_activations = get_activations('masks/mask.nii', event, data, metadata, use_nosmooth);
end

% For each voxel, build a mask that's a sphere around it
% Make sure to only include voxels from the brian
%
for i=1:length(x)

    % Build spherical mask
    %
    sphere_mask = create_spherical_mask_helper(mask, x(i), y(i), z(i), r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);

    % Get sphere center coordinates in MNI space
    %
    mni = cor2mni([x(i) y(i) z(i)], Vmask.mat);

    if use_pregenerated_activations
        sphere_activations = get_activations_submask(sphere_mask, whole_brain_activations); % fast
    else
        sphere_activations = load_activations(sphere_mask, event, data, metadata, use_nosmooth); % old fashioned SUPER SLOW
    end

    searchlight_idx = searchlight_idx + 1;
    Searchlight(searchlight_idx).activations = sphere_activations;
    Searchlight(searchlight_idx).name = ['sphere_', sprintf('%d_%d_%d', mni), '_', event(1)];
    Searchlight(searchlight_idx).center = [x(i) y(i) z(i)];
    Searchlight(searchlight_idx).radius = r;
    Searchlight(searchlight_idx).center_mni = mni;
end


disp('Computed searchlights.');
toc
