% Idea: RDMs of ROIs vs. models (e.g. posterior) on a run-by-run basis;
% correlate with behavior 
% same as kl_structure_learning effect but multivoxel pattern
%

close all;
clear all;

utils;

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

%% Get the neural RDMs
%
%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures', 0.001, '+', 0.001, 20, 1, 1.814);
Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, 'rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', 0.001, '+', 0.001, 20, 1, 1.814);

%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, 'rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', 0.001, '+', 0.001, 20, 1, 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, 'rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', 0.001, '+');
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures', 0.001, '+');
%Neural = rdms_get_glm_and_searchlight_rois(data, metadata, which_trials);
%Neural_controls = rdms_get_anatomical_rois(data, metadata, which_trials, false, false);
%Neural = [Neural, Neural_controls];
%showRDMs(Neural, 1);
Neural = Neural(numel(Neural)/2+1:end); % cut the trail_onset bs


%% Get the model RDMs
%
Model = rdms_get_model(data, metadata, which_trials);
%showRDMs(Model, 2);

control_model_idxs = [8, 12]; % #KNOB control for time and run
assert(isequal(Model(8).name, 'time'));
assert(isequal(Model(12).name, 'run'));

model_idx = 1;
assert(isequal(Model(1).name, 'posterior'));


%% find how similar representations are to posterior (model_idx = 1) in different ROIs per-run
%
rhos = nan(metadata.N, metadata.runsPerSubject, numel(Neural));

goodSubjects = getGoodSubjects();
subjs = metadata.allSubjects(goodSubjects);

% t1, t2 are trial indices corresponding to the pairs of trials in each
% cell of the RDM (for a single subject)
% used to quickly generate RDM masks and extract sub-RDMs.
% Assumes the same trials were used from all subjects (so it takes subj 1
% for convenience)
%
which_trials_per_subj = which_trials & strcmp(data.participant, subjs{1});
[t1, t2] = meshgrid(find(which_trials_per_subj), find(which_trials_per_subj));

% for each ROI,
% for each run of each subject, get the average RDM
%
for run = 1:metadata.runsPerSubject
    t1_mask = data.runId(t1) == run & data.trialId(t1);
    t2_mask = data.runId(t2) == run & data.trialId(t2);
    run_mask = t1_mask & t2_mask & t1 > t2;
    run_mask_control = repmat(run_mask, 1, 1, numel(control_model_idxs));
    
    for subj = 1:metadata.N
        fprintf('subject %d, run %d\n', subj, run);
        
        for neural_idx = 1:numel(Neural)
            % For a given ROI,
            % For each run for each subject,
            % use a binary mask that says which pairs of trials from the RDM
            % to look at.
            %
            neural_RDM = Neural(neural_idx).RDMs(:,:,subj);
            model_RDM = Model(model_idx).RDMs(:,:,subj);
            control_RDMs = nan(size(neural_RDM, 1), size(model_RDM, 2), numel(control_model_idxs));
            for i = 1:numel(control_model_idxs) 
                control_RDMs(:,:,i) = Model(control_model_idxs(i)).RDMs(:,:,subj);
            end
            
            neural_subRDM = neural_RDM(run_mask);
            model_subRDM = model_RDM(run_mask);
            control_subRDMs = control_RDMs(run_mask_control);
            
            control_subRDMs = reshape(control_subRDMs, size(neural_subRDM, 1), numel(control_model_idxs));
            %rho = partialcorr(neural_subRDM, model_subRDM, control_subRDMs, 'type', 'Spearman');
            rho = corr(neural_subRDM, model_subRDM, 'type', 'Spearman');
            assert(~isnan(rho));

            rhos(subj, run, neural_idx) = rho; 
        end
    end
end

%% Get the test choice likelihoods and plot the correlation
%
test_log_liks = get_test_behavior();

correlate_neural_and_behavior(rhos, {Neural.name}, test_log_liks, 'Neural ~ posterior, correlated w/ test log likelihood: t-test');