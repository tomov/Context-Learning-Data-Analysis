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

[inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs, activations, predict_what, z_score, data, metadata);
