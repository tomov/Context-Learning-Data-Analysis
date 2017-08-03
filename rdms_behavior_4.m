% Idea: compute neural RDMs on the run-level (rather than the trial level)
% from ROIs from the searchlight -> see how similar runs look based on ROIs
% that encode the posterior -> see how similar the encoded posterior might
% be on different runs
% then, compare test trial responses on pairs of runs -- see how
% similar they are (based on x1c1, x1c3, x3c1, x3c3)
% if two runs seem to encode similar posteriors => the subject responses on
% the test trials of these runs should probably be similar too
% => ...
%

close all;
clear all;

% load('results/rdms_behavior.mat'); % <--- alternative to recomputing all
% the shit from scratch

utils;

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

%% get the neural rdms
%
Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, 'rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', 0.001, '+', 0.001, 20, 3, 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, 'rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', 0.001, '+');

%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_weights', 0.001, '+', 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_weights', 0.001, '+');

%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures', 0.001, '+', 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures', 0.001, '+');

%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_weights - KL_structures', 0.001, '+', 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_weights - KL_structures', 0.001, '+');

%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures - KL_weights', 0.001, '+', 1.814);
%Neural = rdms_get_rois_from_contrast(data, metadata, which_trials, context_expt(), 154, 'KL_structures - KL_weights', 0.001, '+');

%Neural = rdms_get_glm_and_searchlight_rois(data, metadata, which_trials);
%Neural = rdms_get_anatomical_rois(data, metadata, which_trials);
%Neural = [Neural, Neural_controls];
Neural = Neural(numel(Neural)/2+1:end); % cut the trial_onset bs
%showRDMs(Neural, 1);


%% find how similar representations are in each ROI at the end of training
%
avg_dists = nan(metadata.N, metadata.runsPerSubject, numel(Neural));

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

% First, compute RDMs on the run level (i.e. comparing pairs of runs) based
% on the test trial responses
%
test_RDMs = nan(metadata.runsPerSubject, metadata.runsPerSubject, metadata.N);

for subj = 1:metadata.N
    fprintf('subject %d\n', subj);
    
    test_responses = nan(metadata.runsPerSubject, metadata.testTrialsPerRun);

    for run1 = 1:metadata.runsPerSubject
        fprintf('  run %d\n', run1);
        % Get average behavioral distance for (subject, run)
        % TODO account for timeouts
        %
        which_test_trials = data.which_rows & ~data.isTrain & data.runId == run1 & strcmp(data.participant, subjs{subj});
        assert(sum(which_test_trials) <= metadata.testTrialsPerRun);
        x1c1 = strcmp(data.response.keys(which_test_trials & data.cueId == 0 & data.contextId == 0), 'left');
        x1c3 = strcmp(data.response.keys(which_test_trials & data.cueId == 0 & data.contextId == 2), 'left');
        x3c1 = strcmp(data.response.keys(which_test_trials & data.cueId == 2 & data.contextId == 0), 'left');
        x3c3 = strcmp(data.response.keys(which_test_trials & data.cueId == 2 & data.contextId == 2), 'left');
        
        test_responses(run1, :) = [x1c1 x1c3 x3c1 x3c3];
    end
    
    test_RDMs(:,:,subj) = squareRDMs(pdist(test_responses, 'hamming')); % what percentage of coordinates differ
end

test_RDM = mean(test_RDMs, 3);
%showRDMs(test_RDMs, 1);

% Then, compute neural RDMs on the run level (i.e. pairs of runs) based on
% training trial activations for different ROIs
%
clear Neural_run;
for i=1:numel(Neural)
    Neural_run(i).RDMs = zeros(metadata.runsPerSubject, metadata.runsPerSubject, metadata.N);
    Neural_run(i).name = Neural(i).name;
    Neural_run(i).color = Neural(i).color;
end

for run1 = 1:metadata.runsPerSubject
    t1_mask = data.runId(t1) == run1 & data.trialId(t1) > 0;

    for run2 = 1:run1 - 1
        t2_mask = data.runId(t2) == run2 & data.trialId(t2) > 0;
        run_mask = t1_mask & t2_mask;
        %assert(sum(run_mask(:)) == metadata.trainingTrialsPerRun ^ 2 / 4);

        for subj = 1:metadata.N
            fprintf('subject %d, runs %d,%d\n', subj, run1, run2);

            for neural_idx = 1:numel(Neural)
                RDM = Neural(neural_idx).RDMs(:,:,subj);
                sub_RDM = RDM(run_mask);
                
                Neural_run(neural_idx).RDMs(run1, run2, subj) = mean(sub_RDM);
                Neural_run(neural_idx).RDMs(run2, run1, subj) = Neural_run(neural_idx).RDMs(run1, run2, subj);
            end
        end
    end
end

for i=1:numel(Neural)
    Neural_run(i).RDM = mean(Neural_run(i).RDMs, 3);
end
%showRDMs(Neural_run, 2);



% Correlate test choices and RDMs
%
[r1, r2] = meshgrid(1:metadata.runsPerSubject, 1:metadata.runsPerSubject);
upper = r1 > r2;

neural = nan(metadata.N, sum(upper(:)), numel(Neural_run));
behavioral = nan(metadata.N, sum(upper(:)));

for subj = 1:metadata.N
    RDM = test_RDMs(:,:,subj);
    behavioral(subj,:) = RDM(upper);
    
    for neural_idx = 1:numel(Neural_run)
        RDM = Neural_run(neural_idx).RDMs(:,:,subj);
        neural(subj,:,neural_idx) = RDM(upper);
    end
end

[means, sems, ps] = correlate_neural_and_behavior(neural, {Neural_run.name}, behavioral, 'Training neural ~ test choices', 'Spearman');