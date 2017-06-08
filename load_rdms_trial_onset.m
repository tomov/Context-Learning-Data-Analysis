function [subjectRDMs, avgSubjectRDM, concatSubjectRDMs] = load_rdms_trial_onset(mask, distance_measure, data, metadata)

% Load all subject RDMs for trial_onset events for a given mask and
% distance measure. Load from disk if already computed, otherwise compute
% and save on disk.
%

[~, maskname, ~] = fileparts(mask);
rdms_filename = fullfile('rdms', ['rdms_trial_onset_', maskname, '_', distance_measure, '.mat']);

if exist(rdms_filename, 'file') ~= 2
    fprintf('Computing RDMs from disk and saving them to %s\n', rdms_filename);
    
    which_rows = data.which_rows;
    
    % Load the neural data
    %
    betas = load_betas_trial_onset(mask, data, metadata);

    % Compute the single-subject RDMs
    %
    subjectRDMs = nan(metadata.runsPerSubject * metadata.trialsPerRun, metadata.runsPerSubject * metadata.trialsPerRun, metadata.N);
    
    subj_idx = 0;
    for who = metadata.subjects
        subj_idx = subj_idx + 1;
        who = who{1};
        
        which_subj_rows = which_rows & strcmp(data.participant, who);
        subjectRDMs(:,:,subj_idx) = squareRDMs(pdist(betas(which_subj_rows, :), distance_measure));
    end
    
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