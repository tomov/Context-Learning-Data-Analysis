function [inputs, targets, which_rows] = classify_get_inputs_and_targets(runs, trials, subjs, mask, predict_what, z_score, event)

% helper method that gets the input betas and the targets for the
% classifier to train / test on
%
% INPUT:
% runs = which runs to use, e.g. 1:9 or 1:8
% trials = which trials to use to train from each run e.g. 1:19 or 1:24
% subjs = subject indices in the subjects array returned by
%         contextGetSubjectsDirsAndRuns, e.g. getGoodSubjects()
% mask = .nii file name of which mask to use
% predict_what = 'condition', 'response', 'runId', or 'contextId'
% z_score = how to z-score the betas within each run to account for
%           stuff like drift and faciliate cross-run comparison. Options:
%           'z-none' = no z-scoring,
%           'z-run' = z-score all voxels within each run
%           'z-run-voxel' = z-score each voxel separately within each run
% event = 'feedback_onset' or 'trial_onset'
% 
% OUTPUT:
% inputs = [n_observations x n_voxels] input to classifier
% targets = [n_observations x n_classes] targets for classifier
% which_rows = which trials of the subject data do these inputs correspond
%              to

%fprintf('classify_get_inputs_and_targets\n');
%disp(runs)
%disp(trials)
%disp(mask)
%disp(subjs)
%disp(predict_what);

use_tmaps = false; % #KNOB TODO make parameter
use_nosmooth = false; % #KNOB TODO make parameter

assert(ismember(event, {'trial_onset', 'feedback_onset'}));

if use_tmaps
    get_activations = @get_tmaps;
else
    get_activations = @get_betas;
end

[~,Vwhole] = load_mask(fullfile('masks', 'mask.nii')); % for sanity check

% Load the subject behavioral data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
subjects = metadata.allSubjects;

% Load the neural data
%
[~, name] = system('hostname');
if isempty(strfind(name, 'omchil')) % it could be NCF
    activations = get_activations(mask, event, data, metadata, use_nosmooth); % <-- get the mask straight up, or ...
else % certainly not NCF
    whole_brain_activations = get_activations(fullfile('masks', 'mask.nii'), event, data, metadata, use_nosmooth); % get the submask from the whole-brain betas
    [mask, Vmask] = load_mask(mask);
    assert(isequal(Vwhole.mat, Vmask.mat)); % it is imperative that they are in the same coordinate space if we're getting the betas like this!!!
    activations = get_activations_submask(mask, whole_brain_activations);
end

% condition = context role labels
%
condition_labels = containers.Map(metadata.conditions, ...
                        {[1 0 0], [0 1 0], [0 0 1]});
                    
n_observations = length(subjs) * length(runs) * length(trials);
n_voxels = size(activations, 2);

% figure out which rows (trials) we're looking at
%
which_rows = ismember(data.participant, metadata.allSubjects(subjs)) & ...
    ismember(data.runId, runs) & ismember(data.newTrialId, trials);
assert(~any(data.drop(which_rows)), 'we should only be using good subjects here; maybe passed some "bad" subjects?');
assert(sum(which_rows) == n_observations, 'maybe wrong order of parameters, e.g. runs and trials');

if strcmp(predict_what, 'responses')
    % if we're looking at subject responses, ignore trials when the subject
    % failed to respond
    which_rows = which_rows & ~data.timeout;
    n_observations = sum(which_rows);
end

%
% Compute input vectors
%

inputs = activations(which_rows, :); % rows = x = observations, cols = voxels / dependent vars
inputs(isnan(inputs)) = 0; % some of them fall outside the imaging range and are NaNs; don't let those screw up our analysis
assert(size(inputs, 1) == n_observations);
assert(size(inputs, 2) == n_voxels);

%
% Compute target vectors
%

targets = []; % rows = y = observations, cols = indep vars (condition as binary vector)

% target depends on what we're trying to predict
%
switch predict_what
    case 'condition'
        targets = values(condition_labels, data.contextRole(which_rows));
        targets = cell2mat(targets);
        
    case 'response'
        targets = zeros(n_observations, 2);
        targets(sub2ind(size(targets), [1:n_observations]', data.chose_sick(which_rows) + 1)) = 1;
        
    case 'runId'
        targets = zeros(n_observations, 9);
        targets(sub2ind(size(targets), [1:n_observations]', data.runId(which_rows))) = 1;
        
    case 'contextId'
        targets = zeros(n_observations, 3);
        targets(sub2ind(size(targets), [1:n_observations]', data.contextId(which_rows) + 1)) = 1;
        
    case 'contextId-training-only'
        assert(all(data.contextId(which_rows) < 2));
        targets = zeros(n_observations, 2);
        targets(sub2ind(size(targets), [1:n_observations]', data.contextId(which_rows) + 1)) = 1;
        
    otherwise
        assert(false, 'should be one of the above');
end            

assert(size(targets, 1) == size(inputs, 1));

%
% Optionally z-score the inputs
%

% get the attributes of the trials we took the betas from
% note that from now on, we'll use these to generate bitmasks referring to
% rows from the the inputs vector, instead of the ones in data
%
runId = data.runId(which_rows);
newTrialId = data.newTrialId(which_rows);
participant = data.participant(which_rows);

switch z_score
    
    case 'z-none'
        % do nothing
        
    case 'z-run'
        % z-score all voxels within the same run
        %
        for subj = subjs
            for run = runs
                which = strcmp(participant, metadata.allSubjects{subj}) & runId == run;
                assert(sum(which) == numel(trials));
                
                inputs_for_run = inputs(which, :);
                z_scored_inputs_for_run = reshape(zscore(inputs_for_run(:)), size(inputs_for_run));
                inputs(which, :) = z_scored_inputs_for_run;
                assert(abs(mean(z_scored_inputs_for_run(:))) < 1e-10);
            end
        end
        
    case 'z-run-voxel'
        % z-score each voxel within the same run
        %
        for subj = subjs
            for run = runs
                which = strcmp(participant, metadata.allSubjects{subj}) & runId == run;
                assert(sum(which) == numel(trials));
                
                inputs_for_run = inputs(which, :);
                inputs(which, :) = zscore(inputs_for_run, 0, 1);
                assert(max(abs(mean(inputs(which, :), 1))) < 1e-10);
            end
        end
        
    otherwise
        assert(false, 'invalid z_score -- should be one of the above');        
    
end

