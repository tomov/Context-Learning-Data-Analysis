% RDMs
%

% Load the subject behavioral data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

masks = {fullfile('masks', 'mask.nii'), ...
         fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'v1.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
     
distance_measure = 'cosine';     

for mask = masks
    mask = mask{1}

    % Load the RDMs for that mask
    %
    [subjectRDMs, avgSubjectRDM, concatSubjectRDMs] = load_rdms_trial_onset(mask, distance_measure, data, metadata);
    
    % Display the RDMs
    %
    %showRDMs(concatSubjectRDMs, 2);
    showRDMs(avgSubjectRDM, 1);
end

