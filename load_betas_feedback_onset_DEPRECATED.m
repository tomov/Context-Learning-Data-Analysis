function betas = load_betas_feedback_onset(mask, data, metadata)

% Load all the betas corresponding to activations at feedback onset for the
% given mask. Does this for all trials for all runs for all subjects.
% First try to get the cached betas from disk; if not, loads them and then
% saves the cached file for later use.
%
% INPUT:
% mask = path to .nii file with the mask, e.g. 'masks/hippocampus.nii'
% data, metadata = subject behavioral data as output by load_data
%
% OUTPUT:
% betas = [nRows x nVoxels] beta coefficients for the mask for each trial.
%         Each row of this matrix corresponds to a single trial == a single 
%         row in data.
%

% Load the neural data
%
[~, maskname, ~] = fileparts(mask);
betas_filename = fullfile('betas', ['betas_feedback_onset_', maskname, '.mat']);

% Try to load the precomputed betas from the single mat file. If they
% haven't been precomputed yet, load them and save them in a file for
% future use
%
if exist(betas_filename, 'file') ~= 2
    fprintf('Loading betas from disk and saving them to %s\n', betas_filename);
    betas = load_betas(mask, 'feedback_onset', data, metadata);
    save(betas_filename, 'betas');
else
    fprintf('Loading precomputed betas from %s\n', betas_filename);
    load(betas_filename, 'betas'); % crucial to load betas only
end
