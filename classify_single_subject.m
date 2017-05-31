% Script that train a classifier to predict condition based on activity at trial onset
% separately for each subject.
% Based on classify_single_subject.sh which overwhelmed NCF lol
%

method = 'cvglmnet';
masks = {fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
z_score = 'z-none';
predict_what = 'condition';
runs = [1:9];
trials = [1:24];

subjects = getGoodSubjects();

tic

for mask = masks
    mask = mask{1};
    [~, maskname, ~] = fileparts(mask);
    
    fprintf('\n\n----- Mask %s ---\n\n', maskname);
    
    filename = fullfile('classifier', ['single_subject_', method, '_', z_score, '_', predict_what, '_', maskname, '.mat']);
    fprintf('filename = %s\n', filename);
    
    classifiers = {}; % classifier for each subject
    for subj = subjects
        fprintf('  ---- Subject %d ---------\n', subj);
        
        classifier = classify_train(method, runs, trials, subj, mask, predict_what, z_score);
        classifiers{subj} = classifier;
    end
    
    fprintf('Saving classifiers to %s\n', filename);
    save(filename, '-v7.3');
end
    
toc;