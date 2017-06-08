function [subjectRDMs, avgSubjectRDM, concatSubjectRDMs] = get_neural_rdms(mask, distance_measure, regressor_prefix, data, metadata)

% Compute all subject RDMs for trial_onset or feedback_onset events for a given mask and
% distance measure. Load from disk if already computed, otherwise compute
% and save on disk.

% INPUT:
% mask = path to .nii file with the mask, e.g. 'masks/hippocampus.nii'
% distance_measure = distance argument to pdist, e.g. 'correlation' or
%                    'cosine' or 'euclidean'
% regressor_prefix = 'trial_onset' or 'feedback_onset'; assumes regressor
%                    names are of the form {regressor prefix}_{trial id}
% data, metadata = subject behavioral data as output by load_data
%
%

[~, maskname, ~] = fileparts(mask);
rdms_filename = fullfile('rdms', ['rdms_', regressor_prefix, '_', maskname, '_', distance_measure, '.mat']);

if exist(rdms_filename, 'file') ~= 2
    fprintf('Computing RDMs from disk and saving them to %s\n', rdms_filename);

    betas = get_betas(mask, regressor_prefix, data, metadata);
    subjectRDMs = compute_rdms(betas, distance_measure, data, metadata, data.which_rows);

    save(rdms_filename, 'subjectRDMs', '-v7.3');
else
    fprintf('Loading precomputed RDMs from %s\n', rdms_filename);
    load(rdms_filename, 'subjectRDMs'); % crucial to load betas only
end

% Compute the average subject RDM
% display with showRDMs(avgSubjectRDM, 1)
%
avgSubjectRDM = mean(subjectRDMs, 3);

% For displaying with showRDMs(concatSubjectRDMs, 2)
%
concatSubjectRDMs = cat(3, subjectRDMs, avgSubjectRDM);
