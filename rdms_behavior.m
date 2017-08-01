% Idea: the posterior (and learning in general) is pretty flat in the
% second half of training
% => areas which encode the posterior should look similar during the second
% half of training (as opposed to areas who don't / during the first half
% of training)
% => the extent to which these representations are consistent is indicative
% of how well the subject has learned (...or maybe is not learning at
% all??? idk...)
% anyway, see if there's a relationship between the similarity of
% representations in given ROIs during the second half of training, and
% performance on the test trials (the structure learning effect)
%
% => GLM1 has the effect; however could be just a reflection of the structure
% learning effect we saw earlier...

close all;
clear all;

% load('results/rdms_behavior.mat'); % <--- alternative to recomputing all
% the shit from scratch

utils;

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

%% get the neural rdms
%
Neural = rdms_get_glm_and_searchlight_rois(data, metadata, which_trials);
%Neural = rdms_get_anatomical_rois(data, metadata, which_trials);
%Neural = [Neural, Neural_controls];
showRDMs(Neural, 1);

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

% for each ROI,
% for each run of each subject, get the average RDM
%
for run = 1:metadata.runsPerSubject
    t1_mask = data.runId(t1) == run & data.trialId(t1) > 10;
    t2_mask = data.runId(t2) == run & data.trialId(t2) > 10;
    run_mask = t1_mask & t2_mask & t1 > t2;
    
    for subj = 1:metadata.N
        fprintf('subject %d, run %d\n', subj, run);
        
        for neural_idx = 1:numel(Neural)
            % For a given ROI,
            % For each run for each subject,
            % use a binary mask that says which pairs of trials from the RDM
            % to look at.
            %
            RDM = Neural(neural_idx).RDMs(:,:,subj);
            sub_RDM = RDM(run_mask);
            
            % TODO is this legit? Fisher z-transform?
            avg_dists(subj, run, neural_idx) = mean(sub_RDM);
        end
    end
end

%% Get the test choice likelihoods and plot the correlation
%
test_log_liks = get_test_behavior();

correlate_neural_and_behavior(avg_dists, {Neural.name}, test_log_liks, 'Stability during trials 11..20 correlated w/ test log likelihood: t-test');