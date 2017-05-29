function classifier = classify(method, mask, what, training_runs, training_trials, test_runs, test_trials)

% Train and test classifier to predict aspects of the task from the neural
% data i.e. P(regressor | neural activity)
%
% USAGE example:
% classifier = classify('cvglmnet', 'masks/hippocampus.nii', 'condition', 1:9, 1:19, 1:9, 20:20)
%
% INPUT
% method = which classifier method to use. Supported:
%          'cvglmnet' (recommended), 'glmnet', or 'patternnet'
% mask = path to .nii file with which mask to use, e.g.
%        'masks/hippocampus.nii'
% what = what to classify. Supported:
%        'condition', 'response', 'runId', or 'contextId'
% training_runs = which runs to use to train the classifier, e.g. 1:9 or 1:8
% training_trials = which trials to use to train from each run e.g. 1:19 or 1:24
% test_runs = which runs to use to test the classifier, e.g. 1:9 or 9:9
% test_trials = which trials to use to test from each run, e.g. 20:20 or 1:24
%

subjs = getGoodSubjects();
classifier = classify_train(method, training_runs, training_trials, subjs, mask, what, true);
classify_test(method, classifier, test_runs, test_trials, subjs, mask, what, true);