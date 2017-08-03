function Model = rdms_get_model_2(data, metadata, which_rows)

% Compute the RDMs for some more different models
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
%
% OUTPUT:
% Model = struct array of RDMs

disp('Computing model RDMs...');
tic

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

% Posterior: normalized correlation of posterior
%
for which_one = 1:3
    which_structs = logical([0 0 0 0]);
    which_structs(which_one) = 1;
    [posteriorRDMs, avgPosteriorRDM] = compute_rdms(simulated.P(:, which_structs), 'euclidean', data, metadata, which_rows);
    model_idx = model_idx + 1;
    Model(model_idx).RDMs = posteriorRDMs;
    Model(model_idx).RDM = avgPosteriorRDM;
    Model(model_idx).name = ['posterior', num2str(which_one)];
    Model(model_idx).color = [0 1 0];
end


% Time = seconds since start of run: -abs(T1 - T2)
%
trial_onset = cellfun(@str2double, data.actualChoiceOnset);
[timeRDMs, avgTimeRDM] = compute_rdms(trial_onset, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = timeRDMs;
Model(model_idx).RDM = avgTimeRDM;
Model(model_idx).name = 'time';
Model(model_idx).color = [0 1 0];

% run: 1 if same, 0 o/w
%
[runRDMs, avgRunRDM] = compute_rdms(data.runId, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = runRDMs;
Model(model_idx).RDM = avgRunRDM;
Model(model_idx).name = 'run';
Model(model_idx).color = [0 1 0];



disp('Computed model RDMs.');
toc
