function Neural = rdms_get_spheres_from_contrast(data, metadata, which_rows, EXPT, model, contrast, p, direct, alpha, Dis, Num, r)
% Compute the neural RDMs for a bunch of ROIs based on spheres around the peak voxels of clusters from the
% GLM. Input is same as create_sphere_masks_from_contrast
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
% EXPT = experiment structure, e.g. context_expt()
% model = GLM number, e.g. 154
% contrast = contrast, e.g. 'KL_weights - KL_structures'
% p = optional p-value threshold
% direct = sign of the activations; should be one of +, -, or +/-
% alpha = significance level for cluster FWE correction, default 0.05 in bspmview
% Dis = separation for cluster maxima, default 20 in bspmview
% Num = numpeaks for cluster maxima, default 3 in bspmview
% r = sphere radius in voxels (native coordinate space)
%
% OUTPUT:
% Neural = struct array of RDMs

assert(ismember(direct, {'+/-', '+', '-'}));

events = {'trial_onset', 'feedback_onset'};
use_tmaps = false;
use_nosmooth = false;

[filenames, masknames] = create_sphere_masks_from_contrast(EXPT, model, contrast, p, direct, alpha, Dis, Num, r);

for mask_idx = 1:numel(filenames)
    masks(mask_idx).filename = filenames{mask_idx};
    masks(mask_idx).rdm_name = masknames{mask_idx};
end

Neural = rdms_get_neural(masks, events, data, metadata, which_rows, use_tmaps, use_nosmooth);