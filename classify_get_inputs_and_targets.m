function [inputs, targets, which_rows] = classify_get_inputs_and_targets(runs, trials, subjs, mask, predict_what, z_score)

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
% z_score = whether to z-score the betas within each run to account for
%           stuff like drift and faciliate cross-run comparison
% 
% OUTPUT:
% inputs = [n_observations x n_voxels] input to classifier
% targets = [n_observations x n_classes] targets for classifier
% which_rows = which trials of the subject data do these inputs correspond
%              to

fprintf('classify_get_inputs_and_targets\n');
disp(runs)
disp(trials)
disp(mask)
disp(subjs)
disp(predict_what);

% Load the subject behavioral data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
subjects = metadata.allSubjects;

% Load the neural data
%
[~, maskname, ~] = fileparts(mask);
betas_filename = fullfile('betas', ['betas_trial_onset_', maskname, '.mat']);

% Try to load the precomputed betas from the single mat file. If they
% haven't been precomputed yet, load them and save them in a file for
% future use
%
if exist(betas_filename, 'file') ~= 2
    fprintf('Loading betas from disk and saving them to %s\n', betas_filename);
    betas = load_betas_trial_onset(mask, data, metadata);
    save(betas_filename, '-v7.3');
else
    fprintf('Loading precomputed betas from %s\n', betas_filename);
    load(betas_filename);
end

% condition = context role labels
%
condition_labels = containers.Map(metadata.conditions, ...
                        {[1 0 0], [0 1 0], [0 0 1]});
                    
n_observations = length(subjs) * length(runs) * length(trials);
n_voxels = size(betas, 2);

which_rows = ~data.drop & ismember(data.participant, metadata.allSubjects(subjs)) & ...
    ismember(data.runId, runs) & ismember(data.newTrialId, trials);
save('test.mat');
assert(sum(which_rows) == n_observations);

if strcmp(predict_what, 'responses')
    % if we're looking at subject responses, ignore trials when the subject
    % failed to respond
    which_rows = which_rows & ~data.timeout;
    n_observations = sum(which_rows);
end

% Compute input vectors
%
inputs = betas(which_rows, :); % rows = x = observations, cols = voxels / dependent vars
assert(size(inputs, 1) == n_observations);
assert(size(inputs, 2) == n_voxels);

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
    case 'contextId_training_only'
        assert(all(data.contextId(which_rows) < 2));
        targets = zeros(n_observations, 2);
        targets(sub2ind(size(targets), [1:n_observations]', data.contextId(which_rows) + 1)) = 1;
end            

assert(size(targets, 1) == size(inputs, 1));


%{
if z_score
    % TODO
    run_mean_voxel = mean(reshape(inputs(run_input_idxs, :), 1, length(run_input_idxs) * n_voxels));
    run_std_voxel = std(reshape(inputs(run_input_idxs, :), 1, length(run_input_idxs) * n_voxels));
    inputs(run_input_idxs, :) = (inputs(run_input_idxs, :) - run_mean_voxel) / run_std_voxel;
end
%}