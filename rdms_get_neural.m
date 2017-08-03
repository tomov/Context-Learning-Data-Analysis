function Neural = rdms_get_neural(masks, events, data, metadata, which_rows, use_tmaps, use_nosmooth)

% Compute the neural RDMs
% Normalized correlation for neural data for different ROIs
%
% INPUT:
% masks = struct array with the masks that we're using for the RDMs
%     .filename = path to .nii file with the mask
%     .rdm_name = name of the RDM; suffix will be added corresponding to
%                 each event for which the RDM is computed
% events = cell array of event names, e.g. {'trial_onset', 'feedback_onset'}.
%          Each mask gets a separate RDM for each event
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

[~,Vwhole] = load_mask(fullfile('masks', 'mask.nii')); % for sanity check

% for each ROI, take betas at both trial onset and feedback onset
%
for event = events
    event = event{1};
    assert(ismember(event, {'trial_onset', 'feedback_onset'}));

    % Load neural data
    %
    whole_brain_activations = get_activations(fullfile('masks', 'mask.nii'), event, data, metadata, use_nosmooth);
    
    for i = 1:numel(masks)
        [mask, Vmask] = load_mask(masks(i).filename);
        assert(isequal(Vwhole.mat, Vmask.mat)); % it is imperative that they are in the same coordinate space if we're getting the betas like this !!!

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
