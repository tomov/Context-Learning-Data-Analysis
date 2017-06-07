function [classifiers, all_targets, all_outputs] = classify_each_subject_CV(mask, z_score)

% Train a separate classifier for each subject to predict the run condition
% based on neural activation.
% Cross-validate using leave-one-run-out to assess accuracy. Have a separate classifier for each
% fold.
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
classifiers_filename = fullfile('classifier', ['each_subject_CV_', method, '_', z_score, '_', predict_what, '_', maskname, '.mat']);
fprintf('filename = %s\n', classifiers_filename);

if exist(classifiers_filename, 'file') ~= 2

    classifiers = {}; % classifier for each fold for each subject
    all_targets = []; % classifier targets *on left out trials* for all subjects / rows of the data
    all_outputs = []; % classifier outputs *on left out trials* for all subjects / rows of the data
    all_which_rows = logical([]); % which rows were used
    accuracies = {}; % accuracy for each fold for each subject
    for subj = subjects
        fprintf('  ---- Classifying subject %d ---------\n', subj);

        fold = 0;
        for leftout_run = runs
            training_runs = runs(runs ~= leftout_run);
            fold = fold + 1;
            fprintf('         ---- Fold %d -- left out run %d, training runs =\n', fold, leftout_run);
            disp(training_runs);
            
            [classifier] = classify_train(method, training_runs, trials, subj, mask, predict_what, z_score);
            classifiers{fold}{subj} = classifier;
            
            [~, targets, outputs, which_rows] = classify_test(method, classifier, leftout_run, trials, subj, mask, predict_what, z_score);
            all_targets(which_rows, :) = targets;
            all_outputs(which_rows, :) = outputs;
            all_which_rows(which_rows) = 1;
            
            accuracy = classify_get_accuracy(outputs, targets);
            accuracies{fold}{subj} = accuracy;
            fprintf('             Success rate is %.2f%% (on %d left-out trials)\n', accuracy, sum(which_rows));
        end        
    end

    fprintf('Saving classifiers to %s\n', classifiers_filename);
    save(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', 'all_which_rows', 'accuracies', '-v7.3');
else
    fprintf('Loading saved classifiers\n');
    load(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', 'all_which_rows', 'accuracies');
end

