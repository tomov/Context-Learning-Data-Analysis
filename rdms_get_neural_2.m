function Neural = rdms_get_neural_2(masks, data, metadata, which_rows, use_tmaps, use_nosmooth)

% TODO dedupe with rdms_get_neural -- make that one use this one
%

% Compute the neural RDMs
% Normalized correlation for neural data for different ROIs
%
% INPUT:
% masks = struct array with the masks that we're using for the RDMs
%     .filename = path to .nii file with the mask
%     .rdm_name = name of the RDM; suffix will be added corresponding to
%                 each event for which the RDM is computed
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
% use_tmaps = if true, use tmaps; otherwise, use betas for neural data
% use_nosmooth = whether to use the non-smoothed neural data
%
% OUTPUT:
% Neural = struct array of RDMs

neural_idx = 0;

disp('Computing neural RDMs...');
tic

if use_tmaps
    get_activations = @get_tmaps;
else
    get_activations = @get_betas;
end

% for each ROI, take betas at both trial onset and feedback onset
%
for event = {'trial_onset', 'feedback_onset'}
    event = event{1};

    % Load neural data
    %
    whole_brain_activations = get_activations(fullfile('masks', 'mask.nii'), event, data, metadata, use_nosmooth);
    
    for i = 1:numel(masks)
        mask = load_mask(masks(i).filename);
        activations = get_activations_submask(mask, whole_brain_activations);
        [hippocampusRDMs, avgHippocampusRDM] = compute_rdms(activations, 'cosine', data, metadata, which_rows);
        neural_idx = neural_idx + 1;
        Neural(neural_idx).RDMs = hippocampusRDMs;
        Neural(neural_idx).RDM = avgHippocampusRDM;
        Neural(neural_idx).name = [masks(i).rdm_name, '_', event(1)];
        Neural(neural_idx).color = [0 1 0];
    end
end

disp('Computed neural RDMs.');
toc
