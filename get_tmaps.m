function tmaps = get_tmaps(mask, regressor_prefix, data, metadata)

% Wrapper around load_tmaps that caches the results to a file and tries to load from that file
% subsequently.
% TODO dedupe with get_betas
%

assert(ismember(regressor_prefix, {'trial_onset', 'feedback_onset'}));

% Load the neural data
%
[~, maskname, ~] = fileparts(mask);
tmaps_filename = fullfile('tmaps', ['tmaps_', regressor_prefix, '_', maskname, '.mat']);

% Try to load the precomputed tmaps from the single mat file. If they
% haven't been precomputed yet, load them and save them in a file for
% future use
%
if exist(tmaps_filename, 'file') ~= 2
    fprintf('Loading tmaps from disk and saving them to %s\n', tmaps_filename);
    tmaps = load_tmaps(mask, regressor_prefix, data, metadata);
    save(tmaps_filename, 'tmaps', '-v7.3');
else
    fprintf('Loading precomputed tmaps from %s\n', tmaps_filename);
    load(tmaps_filename, 'tmaps'); % crucial to load tmaps only
end
