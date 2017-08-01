function Neural = rdms_get_anatomical_rois(data, metadata, which_rows)
% Compute the neural RDMs for a bunch of anatomically defined ROIs
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
%
% OUTPUT:
% Neural = struct array of RDMs

events = {'trial_onset', 'feedback_onset'};
use_tmaps = false;
use_nosmooth = false;

mask_idx = 0;

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/hippocampus.nii';
masks(mask_idx).rdm_name = 'hippocampus';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/ofc.nii';
masks(mask_idx).rdm_name = 'OFC';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/med_ofc.nii';
masks(mask_idx).rdm_name = 'mOFC';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/vmpfc.nii';
masks(mask_idx).rdm_name = 'vmPFC';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/striatum.nii';
masks(mask_idx).rdm_name = 'Striatum';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/pallidum.nii';
masks(mask_idx).rdm_name = 'Pallidum';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/v1.nii';
masks(mask_idx).rdm_name = 'V1';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/m1.nii';
masks(mask_idx).rdm_name = 'M1';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/s1.nii';
masks(mask_idx).rdm_name = 'S1';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/fusiform.nii';
masks(mask_idx).rdm_name = 'Fusiform';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/angular.nii';
masks(mask_idx).rdm_name = 'Angular';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/mid_front.nii';
masks(mask_idx).rdm_name = 'MidFront';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/dl_sup_front.nii';
masks(mask_idx).rdm_name = 'dlSupFront';

Neural = rdms_get_neural(masks, events, data, metadata, which_rows, use_tmaps, use_nosmooth);
