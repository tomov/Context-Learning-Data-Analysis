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

% Training trials only
%
which_rows = data.which_rows & data.isTrain;

%
% Compute the neural RDMs
%

% Load neural data
%
%whole_brain_trial_onset_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata);

% Normalized correlation for neural data for different ROIs
%
clear Neural;
neural_idx = 0;

hippocampus_mask = load_mask('masks/hippocampus.nii');
hippocampus_betas = get_betas_submask(hippocampus_mask, whole_brain_trial_onset_betas);
[hippocampusRDMs, avgHippocampusRDM] = compute_rdms(hippocampus_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgHippocampusRDM;
Neural(neural_idx).name = 'hippocampus';
Neural(neural_idx).color = [0 1 0];

ofc_mask = load_mask('masks/ofc.nii');
ofc_betas = get_betas_submask(ofc_mask, whole_brain_trial_onset_betas);
[ofcRDMs, avgOfcRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgOfcRDM;
Neural(neural_idx).name = 'OFC';
Neural(neural_idx).color = [0 1 0];

vmpfc_mask = load_mask('masks/vmpfc.nii');
vmpfc_betas = get_betas_submask(vmpfc_mask, whole_brain_trial_onset_betas);
[vmpfcRDMs, avgVmpfcRDM] = compute_rdms(vmpfc_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgVmpfcRDM;
Neural(neural_idx).name = 'vmPFC';
Neural(neural_idx).color = [0 1 0];

striatum_mask = load_mask('masks/striatum.nii');
striatum_betas = get_betas_submask(striatum_mask, whole_brain_trial_onset_betas);
[striatumRDMs, avgStriatumRDM] = compute_rdms(striatum_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgStriatumRDM;
Neural(neural_idx).name = 'Striatum';
Neural(neural_idx).color = [0 1 0];

pallidum_mask = load_mask('masks/pallidum.nii');
pallidum_betas = get_betas_submask(pallidum_mask, whole_brain_trial_onset_betas);
[pallidumRDMs, avgPallidumRDM] = compute_rdms(pallidum_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgPallidumRDM;
Neural(neural_idx).name = 'Pallidum';
Neural(neural_idx).color = [0 1 0];

v1_mask = load_mask('masks/v1.nii');
v1_betas = get_betas_submask(v1_mask, whole_brain_trial_onset_betas);
[v1RDMs, avgV1RDM] = compute_rdms(v1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgV1RDM;
Neural(neural_idx).name = 'V1';
Neural(neural_idx).color = [0 1 0];

m1_mask = load_mask('masks/m1.nii');
m1_betas = get_betas_submask(m1_mask, whole_brain_trial_onset_betas);
[m1RDMs, avgM1RDM] = compute_rdms(m1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgM1RDM;
Neural(neural_idx).name = 'M1';
Neural(neural_idx).color = [0 1 0];

s1_mask = load_mask('masks/s1.nii');
s1_betas = get_betas_submask(s1_mask, whole_brain_trial_onset_betas);
[s1RDMs, avgS1RDM] = compute_rdms(s1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgS1RDM;
Neural(neural_idx).name = 'S1';
Neural(neural_idx).color = [0 1 0];

fusiform_mask = load_mask('masks/fusiform.nii');
fusiform_betas = get_betas_submask(fusiform_mask, whole_brain_trial_onset_betas);
[fusiformRDMs, avgFusiformRDM] = compute_rdms(fusiform_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgFusiformRDM;
Neural(neural_idx).name = 'Fusiform';
Neural(neural_idx).color = [0 1 0];

angular_mask = load_mask('masks/angular.nii');
angular_betas = get_betas_submask(angular_mask, whole_brain_trial_onset_betas);
[angularRDMs, avgAngularRDM] = compute_rdms(angular_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgAngularRDM;
Neural(neural_idx).name = 'Angular';
Neural(neural_idx).color = [0 1 0];

mid_front_mask = load_mask('masks/mid_front.nii');
mid_front_betas = get_betas_submask(mid_front_mask, whole_brain_trial_onset_betas);
[midFrontRDMs, avgMidFrontRDM] = compute_rdms(mid_front_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgMidFrontRDM;
Neural(neural_idx).name = 'Middle Frontal';
Neural(neural_idx).color = [0 1 0];

dl_sup_front_mask = load_mask('masks/dl_sup_front.nii');
dl_sup_front_betas = get_betas_submask(dl_sup_front_mask, whole_brain_trial_onset_betas);
[dlSupFrontRDMs, avgDlSupFrontRDM] = compute_rdms(dl_sup_front_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDM = avgDlSupFrontRDM;
Neural(neural_idx).name = 'Superior Frontal, dorsolateral';
Neural(neural_idx).color = [0 1 0];

% show the neural RDMs
% 
showRDMs(Neural, 1);

%
% Compute the Model RDMs
%

simulated = simulate_subjects(data, metadata, params, which_structures);

clear Model;
model_idx = 0;

% Posterior: normalized correlation of posterior
%
[posteriorRDMs, avgPosteriorRDM] = compute_rdms(simulated.P(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'posterior';
Model(model_idx).color = [0 1 0];

% Posterior-ranking: Spearman correlation of posterior
%
[posteriorRankingRDMs, avgPosteriorRankingRDM] = compute_rdms(simulated.P(:, which_structures), 'spearman', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgPosteriorRankingRDM;
Model(model_idx).name = 'posterior-ranking';
Model(model_idx).color = [0 1 0];

% Log posterior: normalized correlation of log(posterior)
%
logPosterior = log(simulated.P + 0.001); % so we don't take log(0);
[logPosteriorRDMs, avgLogPosteriorRDM] = compute_rdms(logPosterior(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgLogPosteriorRDM;
Model(model_idx).name = 'log-posterior';
Model(model_idx).color = [0 1 0];

% Entropy: -abs(H1 - H2) where H = entropy of posterior
%
entropy = - sum(simulated.P(:, which_structures) .* log(simulated.P(:, which_structures)), 2);
entropy(isnan(entropy)) = 0; % if a posterior is 0, the entropy is 0
[entropyRDMs, avgEntropyRDM] = compute_rdms(entropy, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgEntropyRDM;
Model(model_idx).name = 'entropy';
Model(model_idx).color = [0 1 0];

% MAP = Maximum a posteriori: 1 if same structure, 0 o/w
%
[~, map] = max(simulated.P(:, which_structures), [], 2);
[mapRDMs, avgMapRDM] = compute_rdms(map, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgMapRDM;
Model(model_idx).name = 'MAP';
Model(model_idx).color = [0 1 0];

% p(MAP) = probability of MAP structure: -abs(P1 - P2)
%
pMap = max(simulated.P(:, which_structures), [], 2);
[pMapRDMs, avgPMapRDM] = compute_rdms(pMap, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgPMapRDM;
Model(model_idx).name = 'p(MAP)';
Model(model_idx).color = [0 1 0];

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
model_idx = model_idx + 1;
Model(model_idx).RDM = avgPosteriorMapOnlyRDM;
Model(model_idx).name = 'posterior-MAPonly';
Model(model_idx).color = [0 1 0];


% Time = seconds since start of run: -abs(T1 - T2)
%
trial_onset = cellfun(@str2double, data.actualChoiceOnset);
[timeRDMs, avgTimeRDM] = compute_rdms(trial_onset, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDM = avgTimeRDM;
Model(model_idx).name = 'time';
Model(model_idx).color = [0 1 0];

% show the model RDMs
% 
showRDMs(Model, 2);

%
% Second-order similarity matrix from the RDMs
% compare neural RDMs with the model RDMs
%

% Compare the different RDMs
%
userOptions.RDMcorrelationType= 'Spearman';
userOptions.analysisName = 'blah';
userOptions.rootPath = '~/Downloads/'; % TODO how to turn off saving the figure?
corrMat = pairwiseCorrelateRDMs({Neural, Model}, userOptions, struct('figureNumber', 3,'fileName',[]));
