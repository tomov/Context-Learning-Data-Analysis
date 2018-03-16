function [classifier, inputs, targets, outputs, which_rows, accuracy] = classify_train(method, runs, trials, subjs, mask, predict_what, z_score, event)
% Train classifier to predict stuff based on neural activity at trial onset
% returns a fitObj that you can pass to glmnetPredict
% or a petternnet net that you can use e.g. like net(inputs)
% for params descriptions, see classify_get_inputs_and_targets
%
% method = 'cvglmnet', 'glmnet', 'patternnet'


%rng('shuffle');
rng default;

fprintf('classify_train\n');
disp(method);

[inputs, targets, which_rows] = classify_get_inputs_and_targets(runs, trials, subjs, mask, predict_what, z_score, event);
% for debugging, pass the 
%
%{
[data,metadata,simulated] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, [1 1 0 1 0]);
which_rows = data.which_rows & ismember(data.participant, metadata.allSubjects(subjs)) & ...
    ismember(data.runId, runs) & ismember(data.newTrialId, trials);
condition_labels = containers.Map(metadata.conditions, ...
                        {[1 0 0], [0 1 0], [0 0 1]});
targets = values(condition_labels, data.contextRole(which_rows));
targets = cell2mat(targets);
p = simulated.P(which_rows,:);
inputs = p + rand(size(p));
%}



[~, maskname, ~] = fileparts(mask);
outFilename = fullfile('classifier', ['classify_train_', method, '_', maskname, '_', predict_what, '_', z_score, '_', random_string(), '.mat']);

[classifier, outputs, accuracy] = classify_train_helper(method, inputs, targets, runs, trials, subjs, outFilename);
