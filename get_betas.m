function betas = get_betas(mask, regressor_prefix, data, metadata, use_nosmooth)

% Wrapper around load_betas that caches the results to a file and tries to load from that file
% subsequently.
%

assert(ismember(regressor_prefix, {'trial_onset', 'feedback_onset'}));

suffix = '';
if use_nosmooth
    suffix = '_nosmooth';
end

% Load the neural data
%
[~, maskname, ~] = fileparts(mask);
%betas_filename = fullfile('betas', ['betas_', regressor_prefix, '_', maskname, suffix, '.mat']);
betas_filename = fullfile('/Volumes/fMRI/context/code/betas', ['betas_', regressor_prefix, '_', maskname, suffix, '.mat']);

% Try to load the precomputed betas from the single mat file. If they
% haven't been precomputed yet, load them and save them in a file for
% future use
%
if exist(betas_filename, 'file') ~= 2
    fprintf('Loading betas from disk and saving them to %s\n', betas_filename);
    betas = load_betas(mask, regressor_prefix, data, metadata, use_nosmooth);
    save(betas_filename, 'betas', '-v7.3');
else
    fprintf('Loading precomputed betas from %s\n', betas_filename);
    load(betas_filename, 'betas'); % crucial to load betas only
end
