function [classifiers, all_targets, all_outputs] = classify_each_subject(mask, z_score)

% Train a separate classifier for each subject to predict the run condition
% based on neural activation.
%

method = 'cvglmnet';
%z_score = 'z-run'; -- passed as argument
predict_what = 'condition';
runs = [1:9];
trials = [1:20];

% assumes we're working with all good subjects
%
subjects = getGoodSubjects();

% Run computation for each mask
%
[~, maskname, ~] = fileparts(mask);
fprintf('\n\n----- Mask %s ---\n\n', maskname);

% Load the classifiers, if they have been trained already
% Or train them anew and save them if they haven't
%
classifiers_filename = fullfile('classifier', ['single_subject_', method, '_', z_score, '_', predict_what, '_', maskname, '.mat']);
fprintf('filename = %s\n', classifiers_filename);

if exist(classifiers_filename, 'file') ~= 2

    classifiers = {}; % classifier for each subject
    all_targets = []; % classifier targets for all subjects / rows of the data
    all_outputs = []; % classifier outputs for all subjects / rows of the data
    for subj = subjects
        fprintf('  ---- Classifying subject %d ---------\n', subj);

        [classifier, ~, targets, outputs, which_rows] = classify_train(method, runs, trials, subj, mask, predict_what, z_score);
        all_targets(which_rows, :) = targets;
        all_outputs(which_rows, :) = outputs;
        classifiers{subj} = classifier;
    end

    fprintf('Saving classifiers to %s\n', classifiers_filename);
    save(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', '-v7.3');
else
    fprintf('Loading saved classifiers\n');
    load(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets');
end

