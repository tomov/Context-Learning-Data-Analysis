% Load behavior
%
[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

% Simulate behavior
%
load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
% safeguards
assert(options.isFmriData == false);
assert(options.fixedEffects == 1);
assert(isequal(options.which_structures, [1 1 1 0]));
which_structures = logical(options.which_structures);

simulated = simulate_subjects(data, metadata, params, which_structures);

% Load neural data
%
%whole_brain_trial_onset_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata);

% Generate ROI masks
%
hippocampus_mask = load_mask('masks/hippocampus.nii');
ofc_mask = load_mask('masks/ofc.nii');

% Get neural data for ROIs
%
hippocampus_betas = get_betas_submask(hippocampus_mask, whole_brain_trial_onset_betas);
ofc_betas = get_betas_submask(ofc_mask, whole_brain_trial_onset_betas);

%
% Compute a bunch of RDMs
%

% Training trials only
%
which_rows = data.which_rows & data.isTrain;

% Normalized correlation for neural data for different ROIs
%
[ofcRDMs, avgOfcRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);
[hippocampusRDMs, avgHippocampusRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);

% Posterior: normalized correlation of posterior
%
[posteriorRDMs, avgPosteriorRDM] = compute_rdms(simulated.P(:, which_structures), 'cosine', data, metadata, which_rows);

% Posterior-ranking: Spearman correlation of posterior
%
[sPosteriorRDMs, avgSPosteriorRDM] = compute_rdms(simulated.P(:, which_structures), 'spearman', data, metadata, which_rows);

% Log posterior: normalized correlation of log(posterior)
%
logPosterior = log(simulated.P + 0.001); % so we don't take log(0);
[logPosteriorRDMs, avgLogPosteriorRDM] = compute_rdms(logPosterior(:, which_structures), 'cosine', data, metadata, which_rows);

% Entropy: -abs(H1 - H2) where H = entropy of posterior
%
entropy = - sum(simulated.P(:, which_structures) .* log(simulated.P(:, which_structures)), 2);
entropy(isnan(entropy)) = 0; % if a posterior is 0, the entropy is 0
[entropyRDMs, avgEntropyRDM] = compute_rdms(entropy, 'euclidean', data, metadata, which_rows);

% MAP = Maximum a posteriori: 1 if same structure, 0 o/w
%
[~, map] = max(simulated.P(:, which_structures), [], 2);
[mapRDMs, avgMapRDM] = compute_rdms(map, @(x1, x2) x1 ~= x2, data, metadata, which_rows);

% p(MAP) = probability of MAP structure: -abs(P1 - P2)
%
pMap = max(simulated.P(:, which_structures), [], 2);
[pMapRDMs, avgPMapRDM] = compute_rdms(pMap, 'euclidean', data, metadata, which_rows);

% Posterior MAP only = posterior zerod for all structures except the MAP:
% normalized correlation
%
P = simulated.P(:, which_structures);
[~, map] = max(P, [], 2);
idx = sub2ind(size(P), [1:size(P,1)]', map);
Z = zeros(size(P));
Z(idx) = P(idx);
P = Z;
[posteriorMapOnlyRDMs, avgPosteriorMapOnlyRDM] = compute_rdms(P, 'cosine', data, metadata, which_rows);

% Time = seconds since start of run: -abs(T1 - T2)
%
trial_onset = cellfun(@str2double, data.actualChoiceOnset);
[timeRDMs, avgTimeRDM] = compute_rdms(trial_onset, 'euclidean', data, metadata, which_rows);