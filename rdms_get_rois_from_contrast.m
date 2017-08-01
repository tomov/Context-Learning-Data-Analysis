function Neural = rdms_get_rois_from_contrast(data, metadata, which_rows, EXPT, model, contrast, p, direct)
% Compute the neural RDMs for a bunch of ROIs based on clusters from the
% GLM. Input is same as create_masks_from_contrast
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
% EXPT = experiment structure, e.g. context_expt()
% model = GLM number, e.g. 154
% contrast = contrast, e.g. 'KL_weights - KL_structures'
% p = optional p-value threshold; defaults to 0.001
% direct = optional sign of the activations; should be one of +, -, or +/-;
%          defaults to +/-
%
% OUTPUT:
% Neural = struct array of RDMs

if ~exist('p', 'var')
    p = 0.001;
end
if ~exist('direct', 'var')
    direct = '+/-';
end
assert(ismember(direct, {'+/-', '+', '-'}));

events = {'trial_onset', 'feedback_onset'};
use_tmaps = false;
use_nosmooth = false;

[filenames, masknames] = create_masks_from_contrast(EXPT, model, contrast, p, direct);

for mask_idx = 1:numel(filenames)
    masks(mask_idx).filename = filenames{mask_idx};
    masks(mask_idx).rdm_name = masknames{mask_idx};
end

Neural = rdms_get_neural(masks, events, data, metadata, which_rows, use_tmaps, use_nosmooth);