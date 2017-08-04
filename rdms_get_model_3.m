function [Model, control_model_idxs] = rdms_get_model_3(data, metadata, which_rows)

% Compute the RDMs for some more different models
% WARNING: a bunch of hardcoded stuffs...
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
%
% OUTPUT:
% Model = struct array of RDMs

disp('Computing model RDMs...');
tic

control_model_idxs = [];

%% Simulate behavior
%
load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);
% OVERRIDE -- use params from pilot data
params = [0.1249 2.0064];
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
% safeguards
assert(options.isFmriData == false);
assert(options.fixedEffects == 1);
assert(isequal(options.which_structures, [1 1 1 0]));
which_structures = logical(options.which_structures);

%% Simulate behavior using Kalman filter
%
simulated = simulate_subjects(data, metadata, params, which_structures);

% vector that specifies the structure corresponding to the condition
%
current_structure = nan(length(data.condition), 1);
for i = 1:length(data.condition)
    if strcmp(data.condition{i}, 'irrelevant')
        current_structure(i) = 1;
    elseif strcmp(data.condition{i}, 'modulatory')
        current_structure(i) = 2;
    else
        assert(strcmp(data.condition{i}, 'additive'));
        current_structure(i) = 3;
    end
end

%
%% Create model RDMs
%


model_idx = 0;

%% models at trial onset
%

% Posterior weights: normalized correlation of posterior weights
%
[posteriorRDMs, avgPosteriorRDM] = compute_rdms(simulated.ww_posterior + 0.001, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'ww_posterior';
Model(model_idx).color = [0 1 0];

% Prior weights: normalized correlation of prior weights
%
[priorRDMs, avgPriorRDM] = compute_rdms(simulated.ww_prior + 0.001, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'ww_prior';
Model(model_idx).color = [0 1 0];


% Posterior weights + Sigma: normalized correlation of posterior weights
% WARNING: hardcoded ...
%
save('shit.mat')
Sigma_posterior = reshape(simulated.Sigma_posterior, 100, 5228)';
[posteriorRDMs, avgPosteriorRDM] = compute_rdms([simulated.ww_posterior Sigma_posterior], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'ww_Sigma_posterior';
Model(model_idx).color = [0 1 0];

% Prior weights + Sigma: normalized correlation of prior weights
% WARNING: hardcoded...
%
Sigma_prior = reshape(simulated.Sigma_prior, 100, 5228)';
[priorRDMs, avgPriorRDM] = compute_rdms([simulated.ww_prior Sigma_prior], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'ww_Sigma_prior';
Model(model_idx).color = [0 1 0];


% Time = seconds since start of run: -abs(T1 - T2)
%
trial_onset = cellfun(@str2double, data.actualChoiceOnset);
[timeRDMs, avgTimeRDM] = compute_rdms(trial_onset, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = timeRDMs;
Model(model_idx).RDM = avgTimeRDM;
Model(model_idx).name = 'time';
Model(model_idx).color = [0 1 0];
control_model_idxs = [control_model_idxs, model_idx]; % this one's a control model

% run: 1 if same, 0 o/w
%
[runRDMs, avgRunRDM] = compute_rdms(data.runId, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = runRDMs;
Model(model_idx).RDM = avgRunRDM;
Model(model_idx).name = 'run';
Model(model_idx).color = [0 1 0];
control_model_idxs = [control_model_idxs, model_idx]; % this one's a control model




disp('Computed model RDMs.');
toc


end

