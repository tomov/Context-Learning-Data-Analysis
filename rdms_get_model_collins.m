function [Model, control_model_idxs, params, which_structures] = rdms_get_model_collins(data, metadata, which_rows)

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


%% Simulate behavior using Kalman filter
%
[params, which_structures] = model_params('results/fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat')
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


% Joint Posterior: normalized correlation of full joint posterior
%
P_C = simulated.posteriors_Zc_and_C;
P_C = reshape(P_C, size(P_C,1)*size(P_C,2), size(P_C,3))';
P_S = simulated.posteriors_Zs_and_S;
P_S = reshape(P_S, size(P_S,1)*size(P_S,2), size(P_S,3))';
[posteriorRDMs, avgPosteriorRDM] = compute_rdms([P_C P_S], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'joint_posterior';
Model(model_idx).color = [0 1 0];

% Joint Prior: normalized correlation of full joint prior
%
P_C = simulated.priors_Zc_and_C;
P_C = reshape(P_C, size(P_C,1)*size(P_C,2), size(P_C,3))';
P_S = simulated.priors_Zs_and_S;
P_S = reshape(P_S, size(P_S,1)*size(P_S,2), size(P_S,3))';
[priorRDMs, avgPriorRDM] = compute_rdms([P_C P_S], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'joint_prior';
Model(model_idx).color = [0 1 0];


% Conditional Posterior: normalized correlation of conditional posterior
%
P_C = simulated.posterior_Zc_given_c;
P_S = simulated.posterior_Zs_given_s;
[posteriorRDMs, avgPosteriorRDM] = compute_rdms([P_C P_S], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'cond_posterior';
Model(model_idx).color = [0 1 0];

% Conditional Prior: normalized correlation of conditional prior
%
P_C = simulated.prior_Zc_given_c;
P_S = simulated.prior_Zs_given_s;
[priorRDMs, avgPriorRDM] = compute_rdms([P_C P_S], 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'cond_prior';
Model(model_idx).color = [0 1 0];


% Q Posterior: normalized correlation of Q-values after the update
%
Q = simulated.posteriors_Q;
Q = reshape(Q, size(Q,1)*size(Q,2), size(Q,3))';
[posteriorRDMs, avgPosteriorRDM] = compute_rdms(Q, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'Q_posteosterior';
Model(model_idx).color = [0 1 0];


% Q Prior: normalized correlation of Q-values before the update
%
Q = simulated.priors_Q;
Q = reshape(Q, size(Q,1)*size(Q,2), size(Q,3))';
[priorRDMs, avgPriorRDM] = compute_rdms(Q, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'Q_prior';
Model(model_idx).color = [0 1 0];



% CONTROL RDMs


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


%showRDMs(Model);

disp('Computed model RDMs.');
toc


end

