% Compare RDMs from different ROIs with model RDMs
%

%% Load data
%

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

% Look at training trials only
%
which_rows = data.which_rows & data.isTrain;

%% ------------------- Compute the neural RDMs ------------------
%

% Normalized correlation for neural data for different ROIs
%
clear Neural;
neural_idx = 0;

% for each ROI, take betas at both trial onset and feedback onset
%
for event = {'trial_onset', 'feedback_onset'}
    event = event{1};

    % Load neural data
    %
    whole_brain_betas = get_betas('masks/mask.nii', event, data, metadata);

    hippocampus_mask = load_mask('masks/hippocampus.nii');
    hippocampus_betas = get_betas_submask(hippocampus_mask, whole_brain_betas);
    [hippocampusRDMs, avgHippocampusRDM] = compute_rdms(hippocampus_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = hippocampusRDMs;
    Neural(neural_idx).RDM = avgHippocampusRDM;
    Neural(neural_idx).name = ['hippocampus_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    ofc_mask = load_mask('masks/ofc.nii');
    ofc_betas = get_betas_submask(ofc_mask, whole_brain_betas);
    [ofcRDMs, avgOfcRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = ofcRDMs;
    Neural(neural_idx).RDM = avgOfcRDM;
    Neural(neural_idx).name = ['OFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    med_ofc_mask = load_mask('masks/med_ofc.nii');
    med_ofc_betas = get_betas_submask(med_ofc_mask, whole_brain_betas);
    [medOfcRDMs, avgMedOfcRDM] = compute_rdms(med_ofc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = medOfcRDMs;
    Neural(neural_idx).RDM = avgMedOfcRDM;
    Neural(neural_idx).name = ['mOFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    vmpfc_mask = load_mask('masks/vmpfc.nii');
    vmpfc_betas = get_betas_submask(vmpfc_mask, whole_brain_betas);
    [vmpfcRDMs, avgVmpfcRDM] = compute_rdms(vmpfc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = vmpfcRDMs;
    Neural(neural_idx).RDM = avgVmpfcRDM;
    Neural(neural_idx).name = ['vmPFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    striatum_mask = load_mask('masks/striatum.nii');
    striatum_betas = get_betas_submask(striatum_mask, whole_brain_betas);
    [striatumRDMs, avgStriatumRDM] = compute_rdms(striatum_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = striatumRDMs;
    Neural(neural_idx).RDM = avgStriatumRDM;
    Neural(neural_idx).name = ['Striatum_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    pallidum_mask = load_mask('masks/pallidum.nii');
    pallidum_betas = get_betas_submask(pallidum_mask, whole_brain_betas);
    [pallidumRDMs, avgPallidumRDM] = compute_rdms(pallidum_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = pallidumRDMs;
    Neural(neural_idx).RDM = avgPallidumRDM;
    Neural(neural_idx).name = ['Pallidum_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    v1_mask = load_mask('masks/v1.nii');
    v1_betas = get_betas_submask(v1_mask, whole_brain_betas);
    [v1RDMs, avgV1RDM] = compute_rdms(v1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = v1RDMs;
    Neural(neural_idx).RDM = avgV1RDM;
    Neural(neural_idx).name = ['V1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    m1_mask = load_mask('masks/m1.nii');
    m1_betas = get_betas_submask(m1_mask, whole_brain_betas);
    [m1RDMs, avgM1RDM] = compute_rdms(m1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = m1RDMs;
    Neural(neural_idx).RDM = avgM1RDM;
    Neural(neural_idx).name = ['M1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    s1_mask = load_mask('masks/s1.nii');
    s1_betas = get_betas_submask(s1_mask, whole_brain_betas);
    [s1RDMs, avgS1RDM] = compute_rdms(s1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = s1RDMs;
    Neural(neural_idx).RDM = avgS1RDM;
    Neural(neural_idx).name = ['S1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    fusiform_mask = load_mask('masks/fusiform.nii');
    fusiform_betas = get_betas_submask(fusiform_mask, whole_brain_betas);
    [fusiformRDMs, avgFusiformRDM] = compute_rdms(fusiform_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = fusiformRDMs;
    Neural(neural_idx).RDM = avgFusiformRDM;
    Neural(neural_idx).name = ['Fusiform_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    angular_mask = load_mask('masks/angular.nii');
    angular_betas = get_betas_submask(angular_mask, whole_brain_betas);
    [angularRDMs, avgAngularRDM] = compute_rdms(angular_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = angularRDMs;
    Neural(neural_idx).RDM = avgAngularRDM;
    Neural(neural_idx).name = ['Angular_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    mid_front_mask = load_mask('masks/mid_front.nii');
    mid_front_betas = get_betas_submask(mid_front_mask, whole_brain_betas);
    [midFrontRDMs, avgMidFrontRDM] = compute_rdms(mid_front_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = midFrontRDMs;
    Neural(neural_idx).RDM = avgMidFrontRDM;
    Neural(neural_idx).name = ['MidFront_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    dl_sup_front_mask = load_mask('masks/dl_sup_front.nii');
    dl_sup_front_betas = get_betas_submask(dl_sup_front_mask, whole_brain_betas);
    [dlSupFrontRDMs, avgDlSupFrontRDM] = compute_rdms(dl_sup_front_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = dlSupFrontRDMs;
    Neural(neural_idx).RDM = avgDlSupFrontRDM;
    Neural(neural_idx).name = ['dlSupFront_', event(1)];
    Neural(neural_idx).color = [0 1 0];
end


% show the neural RDMs
% 
showRDMs(Neural, 1);




%% ------------------ Compute the Model RDMs ----------------------
%

simulated = simulate_subjects(data, metadata, params, which_structures);

clear Model;
model_idx = 0;

tic;

%
% at trial onset
%

% Posterior: normalized correlation of posterior
%
[posteriorRDMs, avgPosteriorRDM] = compute_rdms(simulated.P(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRDMs;
Model(model_idx).RDM = avgPosteriorRDM;
Model(model_idx).name = 'posterior';
Model(model_idx).color = [0 1 0];

% Posterior-ranking: Spearman correlation of posterior
%
[posteriorRankingRDMs, avgPosteriorRankingRDM] = compute_rdms(simulated.P(:, which_structures), 'spearman', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = posteriorRankingRDMs;
Model(model_idx).RDM = avgPosteriorRankingRDM;
Model(model_idx).name = 'postRanking';
Model(model_idx).color = [0 1 0];

% Log posterior: normalized correlation of log(posterior)
%
logPosterior = log(simulated.P + 0.001); % so we don't take log(0);
[logPosteriorRDMs, avgLogPosteriorRDM] = compute_rdms(logPosterior(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logPosteriorRDMs;
Model(model_idx).RDM = avgLogPosteriorRDM;
Model(model_idx).name = 'logPosterior';
Model(model_idx).color = [0 1 0];

% Entropy: -abs(H1 - H2) where H = entropy of posterior
%
entropy = - sum(simulated.P(:, which_structures) .* log(simulated.P(:, which_structures)), 2);
entropy(isnan(entropy)) = 0; % if a posterior is 0, the entropy is 0
[entropyRDMs, avgEntropyRDM] = compute_rdms(entropy, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = entropyRDMs;
Model(model_idx).RDM = avgEntropyRDM;
Model(model_idx).name = 'entropy';
Model(model_idx).color = [0 1 0];

% MAP = Maximum a posteriori: 1 if same structure, 0 o/w
%
[~, map] = max(simulated.P(:, which_structures), [], 2);
[mapRDMs, avgMapRDM] = compute_rdms(map, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = mapRDMs;
Model(model_idx).RDM = avgMapRDM;
Model(model_idx).name = 'MAP';
Model(model_idx).color = [0 1 0];

% p(MAP) = probability of MAP structure: -abs(P1 - P2)
%
pMap = max(simulated.P(:, which_structures), [], 2);
[pMapRDMs, avgPMapRDM] = compute_rdms(pMap, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pMapRDMs;
Model(model_idx).RDM = avgPMapRDM;
Model(model_idx).name = 'pMAP';
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
Model(model_idx).RDMs = posteriorMapOnlyRDMs;
Model(model_idx).RDM = avgPosteriorMapOnlyRDM;
Model(model_idx).name = 'posteriorMAPonly';
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

% food image: 1 if same, 0 o/w
%
food = data.runId * 10 + data.cueId + 1;
[foodRDMs, avgFoodRDM] = compute_rdms(food, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = foodRDMs;
Model(model_idx).RDM = avgFoodRDM;
Model(model_idx).name = 'food';
Model(model_idx).color = [0 1 0];

% restaurant name: 1 if same, 0 o/w
%
restaurant = data.runId * 10 + data.contextId + 1;
[restaurantRDMs, avgRestaurantRDM] = compute_rdms(restaurant, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = restaurantRDMs;
Model(model_idx).RDM = avgRestaurantRDM;
Model(model_idx).name = 'restaurant';
Model(model_idx).color = [0 1 0];

% stimulus = food + restaurant: 1 if same, 0 o/w
%
stimulus = data.runId * 100 + (data.cueId + 1) * 10 + data.contextId + 1;
[stimulusRDMs, avgStimulusRDM] = compute_rdms(stimulus, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = stimulusRDMs;
Model(model_idx).RDM = avgStimulusRDM;
Model(model_idx).name = 'stimulus';
Model(model_idx).color = [0 1 0];

% run: 1 if same, 0 o/w
%
[runRDMs, avgRunRDM] = compute_rdms(data.runId, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = runRDMs;
Model(model_idx).RDM = avgRunRDM;
Model(model_idx).name = 'run';
Model(model_idx).color = [0 1 0];

% condition: 1 if same, 0 o/w
%
m = containers.Map(metadata.conditions, {1, 2, 3});
cond = cellfun(@(x) m(x), data.condition);
[condRDMs, avgCondRDM] = compute_rdms(cond, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = condRDMs;
Model(model_idx).RDM = avgCondRDM;
Model(model_idx).name = 'condition';
Model(model_idx).color = [0 1 0];

% value
%
[valRDMs, avgValRDM] = compute_rdms(simulated.values, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = valRDMs;
Model(model_idx).RDM = avgValRDM;
Model(model_idx).name = 'value';
Model(model_idx).color = [0 1 0];

% log value^2
%
logVal2 = log((simulated.values + 0.001) .^ 2); % so we don't take log(0) or negative
[logVal2RDMs, avgLog2ValRDM] = compute_rdms(logVal2, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logVal2RDMs;
Model(model_idx).RDM = avgLog2ValRDM;
Model(model_idx).name = 'logValueSquared';
Model(model_idx).color = [0 1 0];

% values
%
values = simulated.valuess(:, which_structures);
values = values + rand(size(values)) * 0.001; % so cosine metric works
[valsRDMs, avgValsRDM] = compute_rdms(values, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = valsRDMs;
Model(model_idx).RDM = avgValsRDM;
Model(model_idx).name = 'values';
Model(model_idx).color = [0 1 0];

% log values^2
%
logVals2 = log((simulated.valuess + 0.001) .^ 2); % so we don't take log(0) or negative
[logVals2RDMs, avgLogVals2RDM] = compute_rdms(logVals2, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logVals2RDMs;
Model(model_idx).RDM = avgLogVals2RDM;
Model(model_idx).name = 'logValuesSquared';
Model(model_idx).color = [0 1 0];

% weights before the update
%
ww_before = [simulated.ww1_before simulated.ww2_before simulated.ww3_before];
ww_before = ww_before + rand(size(ww_before)) * 0.001; % so cosine metric works
[weightsBeforeRDMs, avgWeightsBeforeRDM] = compute_rdms(ww_before, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsBeforeRDMs;
Model(model_idx).RDM = avgWeightsBeforeRDM;
Model(model_idx).name = 'weightsBefore';
Model(model_idx).color = [0 1 0];



%
% at feedback onset
%


% D_KL: -abs(D1 - D2) 
%
[surpriseRDMs, avgSurpriseRDM] = compute_rdms(simulated.surprise, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = surpriseRDMs;
Model(model_idx).RDM = avgSurpriseRDM;
Model(model_idx).name = 'KL';
Model(model_idx).color = [0 1 0];

% log(D_KL)
%
logSurprise = log(simulated.surprise + 0.001); % so we don't take log(0)
[logSurpriseRDMs, avgLogSurpriseRDM] = compute_rdms(logSurprise, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logSurpriseRDMs;
Model(model_idx).RDM = avgLogSurpriseRDM;
Model(model_idx).name = 'logKL';
Model(model_idx).color = [0 1 0];

% PE
%
PE = data.outcome - simulated.values;
[peRDMs, avgPeRDM] = compute_rdms(PE, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = peRDMs;
Model(model_idx).RDM = avgPeRDM;
Model(model_idx).name = 'PE';
Model(model_idx).color = [0 1 0];

% PEs
%
PEs = data.outcome - simulated.valuess;
PEs = PEs + rand(size(PEs)) * 0.001; % so cosine metric works
[pesRDMs, avgPesRDM] = compute_rdms(PEs, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pesRDMs;
Model(model_idx).RDM = avgPesRDM;
Model(model_idx).name = 'PEs';
Model(model_idx).color = [0 1 0];

% PE^2
%
PE2 = (data.outcome - simulated.values) .^ 2;
[pe2RDMs, avgPe2RDM] = compute_rdms(PE2, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pe2RDMs;
Model(model_idx).RDM = avgPe2RDM;
Model(model_idx).name = 'PEsquared';
Model(model_idx).color = [0 1 0];

% log(PE^2)
%
logPE2 = log((data.outcome - simulated.values) .^ 2 + 0.001); % so we don't take log(0)
[logPe2RDMs, avgLogPe2RDM] = compute_rdms(logPE2, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logPe2RDMs;
Model(model_idx).RDM = avgLogPe2RDM;
Model(model_idx).name = 'logPEsquared';
Model(model_idx).color = [0 1 0];

% PEs^2
%
PEs2 = (data.outcome - simulated.valuess) .^ 2;
PEs2 = PEs2 + rand(size(PEs2)); % so cosine metric works
[pes2RDMs, avgPes2RDM] = compute_rdms(PEs2, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pes2RDMs;
Model(model_idx).RDM = avgPes2RDM;
Model(model_idx).name = 'PEsSquared';
Model(model_idx).color = [0 1 0];

% log(PEs^2)
%
logPEs2 = log((data.outcome - simulated.valuess) .^ 2 + 0.001); % so we don't take log(0)
[logPes2RDMs, avgLogPes2RDM] = compute_rdms(logPEs2, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logPes2RDMs;
Model(model_idx).RDM = avgLogPes2RDM;
Model(model_idx).name = 'logPEsSquared';
Model(model_idx).color = [0 1 0];

% 1/lambdas = reliability = precision
%
reliability = 1 ./ simulated.lambdas(:, which_structures);
[reliabilityRDMs, avgReliabilityRDM] = compute_rdms(reliability, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = reliabilityRDMs;
Model(model_idx).RDM = avgReliabilityRDM;
Model(model_idx).name = 'reliability';
Model(model_idx).color = [0 1 0];

% likelihoods
%
[likRDMs, avgLikRDM] = compute_rdms(simulated.likelihoods(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = likRDMs;
Model(model_idx).RDM = avgLikRDM;
Model(model_idx).name = 'likelihood';
Model(model_idx).color = [0 1 0];

% log likelihoods
%
logLik = log(simulated.likelihoods(:, which_structures) + 0.001); % so we don't take log(0)
[logLikRDMs, avgLogLikRDM] = compute_rdms(logLik, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logLikRDMs;
Model(model_idx).RDM = avgLogLikRDM;
Model(model_idx).name = 'logLikelihood';
Model(model_idx).color = [0 1 0];

% Outcome
%
[outcomeRDMs, avgOutcomeRDM] = compute_rdms(data.outcome, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = outcomeRDMs;
Model(model_idx).RDM = avgOutcomeRDM;
Model(model_idx).name = 'outcome';
Model(model_idx).color = [0 1 0];

% Correct / wrong
%
[corrRDMs, avgCorrRDM] = compute_rdms(data.response.corr, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = corrRDMs;
Model(model_idx).RDM = avgCorrRDM;
Model(model_idx).name = 'correct';
Model(model_idx).color = [0 1 0];

% RT
%
RT = data.response.rt;
RT(data.timeout) = 3;
[rtRDMs, avgRtRDM] = compute_rdms(RT, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = rtRDMs;
Model(model_idx).RDM = avgRtRDM;
Model(model_idx).name = 'RT';
Model(model_idx).color = [0 1 0];

% Response
%
resp = double(strcmp(data.response.keys, 'left'));
resp(data.timeout) = 2;
[respRDMs, avgRespRDM] = compute_rdms(resp, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = respRDMs;
Model(model_idx).RDM = avgRespRDM;
Model(model_idx).name = 'Response';
Model(model_idx).color = [0 1 0];

% New value
%
[newValRDMs, avgNewValRDM] = compute_rdms(simulated.new_values, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = newValRDMs;
Model(model_idx).RDM = avgNewValRDM;
Model(model_idx).name = 'newValue';
Model(model_idx).color = [0 1 0];

% Log New value
%
logNewVal2 = log((simulated.new_values + 0.001) .^ 2); % so we don't take log(0) or negative
[logNewVal2RDMs, avgLogNewVal2RDM] = compute_rdms(logNewVal2, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logNewVal2RDMs;
Model(model_idx).RDM = avgLogNewVal2RDM;
Model(model_idx).name = 'logNewValueSquared';
Model(model_idx).color = [0 1 0];

% New values -- normalized correlation (cosine) is not a good measure here
%
[newValsRDMs, avgNewValsRDM] = compute_rdms(simulated.new_valuess(:, which_structures), 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = newValsRDMs;
Model(model_idx).RDM = avgNewValsRDM;
Model(model_idx).name = 'newValues';
Model(model_idx).color = [0 1 0];

% value PE = New value - old value
%
valPE = simulated.new_values - simulated.values(1:end-4);
[valPERDMs, avgValPERDM] = compute_rdms(valPE, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = valPERDMs;
Model(model_idx).RDM = avgValPERDM;
Model(model_idx).name = 'valPE';
Model(model_idx).color = [0 1 0];

% value PEs = New values - old values
% normalized correlation (cosine) is not a good measure here
%
valPEs = simulated.new_valuess(:, which_structures) - simulated.valuess(1:end-4, which_structures);
[valPEsRDMs, avgValPEsRDM] = compute_rdms(valPEs, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = valPEsRDMs;
Model(model_idx).RDM = avgValPEsRDM;
Model(model_idx).name = 'valPEs';
Model(model_idx).color = [0 1 0];


% Prior: normalized correlation of prior
%
[priorRDMs, avgPriorRDM] = compute_rdms(simulated.Q(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRDMs;
Model(model_idx).RDM = avgPriorRDM;
Model(model_idx).name = 'prior';
Model(model_idx).color = [0 1 0];

% Prior-ranking: Spearman correlation of prior
%
Qrand = simulated.Q(:, which_structures);
Qrand = Qrand + rand(size(Qrand)) * 0.000001; % to break the ties at trial 1
[priorRankingRDMs, avgPriorRankingRDM] = compute_rdms(Qrand, 'spearman', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorRankingRDMs;
Model(model_idx).RDM = avgPriorRankingRDM;
Model(model_idx).name = 'priorRanking';
Model(model_idx).color = [0 1 0];

% Log prior: normalized correlation of log(prior)
%
logPrior = log(simulated.Q + 0.001); % so we don't take log(0);
[logPriorRDMs, avgLogPriorRDM] = compute_rdms(logPrior(:, which_structures), 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logPriorRDMs;
Model(model_idx).RDM = avgLogPriorRDM;
Model(model_idx).name = 'logPrior';
Model(model_idx).color = [0 1 0];

% Entropy: -abs(H1 - H2) where H = entropy of prior
%
qEntropy = - sum(simulated.Q(:, which_structures) .* log(simulated.Q(:, which_structures)), 2);
qEntropy(isnan(qEntropy)) = 0; % if a prior is 0, the qEntropy is 0
[qEntropyRDMs, avgQEntropyRDM] = compute_rdms(qEntropy, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = qEntropyRDMs;
Model(model_idx).RDM = avgQEntropyRDM;
Model(model_idx).name = 'qEntropy';
Model(model_idx).color = [0 1 0];

% MAQ = Maximum a priori: 1 if same structure, 0 o/w
%
[~, maq] = max(simulated.Q(:, which_structures), [], 2);
[maqRDMs, avgMaqRDM] = compute_rdms(maq, @(x1, x2) x1 ~= x2, data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = maqRDMs;
Model(model_idx).RDM = avgMaqRDM;
Model(model_idx).name = 'MAQ';
Model(model_idx).color = [0 1 0];

% p(MAQ) = probability of MAQ structure: -abs(P1 - P2)
%
pMaq = max(simulated.Q(:, which_structures), [], 2);
[pMaqRDMs, avgPMaqRDM] = compute_rdms(pMaq, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pMaqRDMs;
Model(model_idx).RDM = avgPMaqRDM;
Model(model_idx).name = 'pMAQ';
Model(model_idx).color = [0 1 0];

% Prior MAQ only = prior zerod for all structures except the MAQ:
% normalized correlation
%
Q = simulated.Q(:, which_structures);
[~, maq] = max(Q, [], 2);
idx = sub2ind(size(Q), [1:size(Q,1)]', maq);
Z = zeros(size(Q));
Z(idx) = Q(idx);
Q = Z;
[priorMaqOnlyRDMs, avgPriorMaqOnlyRDM] = compute_rdms(Q, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = priorMaqOnlyRDMs;
Model(model_idx).RDM = avgPriorMaqOnlyRDM;
Model(model_idx).name = 'priorMAQonly';
Model(model_idx).color = [0 1 0];

% weights after the update
%
ww_after = [simulated.ww1_after simulated.ww2_after simulated.ww3_after];
ww_after = ww_after + rand(size(ww_after)) * 0.001; % so cosine metric works
[weightsAfterRDMs, avgWeightsAfterRDM] = compute_rdms(ww_after, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsAfterRDMs;
Model(model_idx).RDM = avgWeightsAfterRDM;
Model(model_idx).name = 'weightsAfter';
Model(model_idx).color = [0 1 0];


% show the model RDMs
% 
showRDMs(Model, 2);

toc;













%
% Second-order similarity matrix from the RDMs
% compare neural RDMs with the model RDMs
%



%% Compare the average RDMs only -- DON'T RUN THIS...
% this is LAME -- look at the one below
%
userOptions.RDMcorrelationType= 'Spearman';
userOptions.analysisName = 'blah';
userOptions.rootPath = '~/Downloads/'; % TODO how to turn off saving the figure?
corrMat = pairwiseCorrelateRDMs({Neural, Model}, userOptions, struct('figureNumber', 3,'fileName',[]));















%% Proper RDM comparison
% compare each neural and model RDMs for subject separately using
% Spearman's rank coefficient, then find which ones are significant
%
tic;

% KNOBS / analysis params
same_run_only = false; % #KNOB compute correlation for each run separately i.e. no comparisons of representations across runs
all_vs_all = false; % #KNOB compare neural vs. neural and model vs. model representations too
control_for_time = true; % #KNOB do partial correlation controlling for the time model
control_for_run = true; % #KNOB 
linear_mixed_effects = true; % KNOB -- do the linear mixed effects (LME) model

lme_neural_idxs = [1 2 5 11 12 14 15 18 24 25]; % which ROIs to consider for LME
lme_model_idxs = [1 3 18 39 41 46]; % which models to consider to LME

assert(~linear_mixed_effects || ~all_vs_all); % can't do both

% Models to control for
%
control_models = [];
if control_for_time
    control_models = [control_models, 8];
    assert(isequal(Model(8).name, 'time'));
end
if control_for_run
    control_models = [control_models, 12];
    assert(isequal(Model(12).name, 'run'));
end

% upper right triangle, excluding diagonal, for all runs
%
cross_run_trig = logical(triu(ones(metadata.runsPerSubject * metadata.trainingTrialsPerRun), 1));
cross_run_trig_control = repmat(cross_run_trig, 1, 1, numel(control_models));

% upper right triangle, excluding diagonal, for a single run
%
same_run_trig = false(metadata.runsPerSubject * metadata.trainingTrialsPerRun);
for run = 1:metadata.runsPerSubject
    s = (run - 1) * metadata.trainingTrialsPerRun + 1;
    e = run * metadata.trainingTrialsPerRun;
    same_run_trig(s:e,s:e) = logical(triu(ones(metadata.trainingTrialsPerRun), 1));
end
same_run_trig_control = repmat(same_run_trig, 1, 1, numel(control_models));


table_Rho = []; % average Spearman's rho for each ROI for each model
table_H = []; % result of hypothesis test for each ROI for each model -- is the correlation significant?
table_P = []; % p-value of hypothesis test for each ROI for each model

all_subject_rhos = nan(numel(Neural), numel(Model), metadata.N);

if all_vs_all
    all = [Neural, Model];
    rows = all;
    cols = all;
else
    rows = Neural; % order here is important for visualizations below and for linear_mixed_effects
    cols = Model;
end

lmes = {};

for row_idx = 1:numel(rows)

    models_subjs_rhos = []; % Spearman rhos: row = model, col = subject
    for col_idx = 1:numel(cols)

        if linear_mixed_effects && ismember(row_idx, lme_neural_idxs) && ismember(col_idx, lme_model_idxs)
            lme_neural = [];
            lme_model = [];
            lme_ids = [];
            lme_controls = [];
        end

        % Compute a Spearman's rank correlation for each subject separately
        %
        subjs_rhos = []; % Spearman rhos: col = subject, one rho per
        for subj = 1:metadata.N
            % get the RDMs to correlate
            %
            row_RDM = rows(row_idx).RDMs(:,:,subj);
            col_RDM = cols(col_idx).RDMs(:,:,subj);
            control_RDMs = nan(size(row_RDM, 1), size(row_RDM, 2), numel(control_models));
            for i = 1:numel(control_models) 
                control_RDMs(:,:,i) = Model(control_models(i)).RDMs(:,:,subj);
            end
            assert(isequal(size(row_RDM), size(col_RDM)));
            assert(isequal(size(cross_run_trig), size(row_RDM)));
            assert(isequal(size(same_run_trig), size(row_RDM)));

            % get the appropriate entries from the RDMs
            %
            if same_run_only
                % Look at RDMs for each run separately => same-run
                % correlations only
                %
                x = row_RDM(same_run_trig);
                y = col_RDM(same_run_trig);
                z = control_RDMs(same_run_trig_control);
            else
                % Look at RDMs for all runs simultaneously => includes
                % cross-run correlations
                %
                x = row_RDM(cross_run_trig);
                y = col_RDM(cross_run_trig);
                z = control_RDMs(cross_run_trig_control);
            end
           
            % compute the correlations
            %
            if numel(control_models) > 0
                z = reshape(z, size(x, 1), numel(control_models));
                subj_rho = partialcorr(x, y, z, 'type', 'Spearman');
                subjs_rhos = [subjs_rhos, subj_rho];
            else
                subj_rho = corr(x, y, 'type', 'Spearman');
                subjs_rhos = [subjs_rhos, subj_rho];
            end
            all_subject_rhos(row_idx, col_idx, subj) = subj_rho;            

            % if we're doing LMEs
            %
            if linear_mixed_effects && ismember(row_idx, lme_neural_idxs) && ismember(col_idx, lme_model_idxs)
                lme_neural = [lme_neural; x];
                lme_model = [lme_model; y];
                lme_ids = [lme_ids; ones(size(x,1),1) * subj];
                lme_controls = [lme_controls; z];
            end
        end
        models_subjs_rhos = [models_subjs_rhos; subjs_rhos]; % for group-level analysis

        % if we're doing LMEs
        %
        if linear_mixed_effects && ismember(row_idx, lme_neural_idxs) && ismember(col_idx, lme_model_idxs)
            assert(isequal(control_models, [8 12])); % time & run
            tbl = array2table([lme_neural, lme_model, lme_ids, lme_controls], 'VariableNames', {'Neural', 'Model', 'Subject', 'Time', 'Run'});
            formula = 'Neural ~ Model + Run + Time + (Model|Subject) + (Run-1|Subject) + (Time|Subject)';
            fprintf('LME for %s vs. %s: %s\n', rows(row_idx).name, cols(col_idx).name, formula);
            lme = fitlme(tbl, formula);
            lmes{row_idx, col_idx} = lme;
        end
    end
    
    % Group-level analysis
    %
    fisher_models_rhos = atanh(models_subjs_rhos);
    [h, ps, ci, stats] = ttest(fisher_models_rhos');

    table_Rho = [table_Rho; mean(fisher_models_rhos')];
    table_H = [table_H; h];
    table_P = [table_P; ps];
end


Rho = array2table(table_Rho, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
H = array2table(table_H, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
P = array2table(table_P, 'RowNames', {rows.name}, 'VariableNames', {cols.name});

toc;
















%% Visualize full analysis
% Show full second-order RSA matrix with correlation coefficients + another one with the p-values
%

tabs = {table_Rho, table_P};
titles = {'Representational similarity match (Spearman rank correlation)', 'P-value (one-sample t-test)'};

for i = 1:numel(tabs)
    figure;
    if i == 2
        t = tabs{i};
        if all_vs_all
            imagesc(log10(t), [-16 0.1]);
        else
            imagesc(log10(t));
        end    
        c = colorbar;
        yt = get(c, 'Ytick');
        set(c, 'YTickLabel', arrayfun(@(x) ['10\^', num2str(x)], yt, 'UniformOutput', false));
    else
        if all_vs_all
            imagesc(tabs{i}, [0 0.4]);
        else
            imagesc(tabs{i});
        end
        %cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
        %colormap(cols); colorbar;
        colorbar;
    end

    xticklabels({cols.name});
    set(gca, 'xtick', 1:numel({cols.name}));
    xtickangle(60);
    xlabel('Neural model');

    yticklabels({rows.name});
    set(gca, 'ytick', 1:numel({rows.name}));
    ylabel('ROI_{event}, t = trial onset, f = feedback onset');
    
    title(titles{i});
end













%% Show the significant positive correlations
%
bonferroni = false; % #KNOB -- whether to use Bonferroni correction for the p-values

if bonferroni
    which = table_Rho > 0 & table_P < 0.05 / numel(table_P); % Bonferroni correction
else
    which = table_Rho > 0 & table_P < 0.01;
end
assert(~all_vs_all); % only works for neural vs. model comparisons

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

plot_idx = 0;
for row_idx = 1:numel(rows)
    col_idxs = find(which(row_idx,:))';
    
    if isempty(col_idxs), continue; end

    % Get fisher z-transformed correlation coefficients
    %
    zs = reshape(all_subject_zs(row_idx, col_idxs, :), [numel(col_idxs) size(all_subject_zs, 3)]);

    % One-sample t-test them against 0
    %
    [h, ps, ci, stats] = ttest(zs');
    assert(numel(h) == numel(col_idxs));

    % Compute SEs
    %
    %[sems, means] = wse(zs'); % WSE's
    means = mean(zs, 2);
    sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % normal SEMs
    assert(numel(sems) == numel(col_idxs));
    assert(numel(means) == numel(col_idxs));

    % Plot
    %
    plot_idx = plot_idx + 1;
    if bonferroni
        subplot(2, 5, plot_idx);
    else
        subplot(3, 7, plot_idx);
    end

    % Plot the bar graphs with error bars
    %
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    % Put the p-values
    %
    for i = 1:numel(col_idxs)
        significance = @(p) repmat('*', 1, floor(-log10(p)) - 1);
        if bonferroni
            stars = significance(ps(i) * numel(table_P)); % Bonferroni correction
        else
            stars = significance(ps(i));
        end
        text(xs(i) - length(stars) * 0.3 / 2, means(i) + sems(i) * 1.2, stars, 'Rotation', 45);
    end

    % Put the ROI / model names and figure title
    %
    set(gca, 'xtick', xs);
    models = {Model(col_idxs).name};
    labels = {};
    for j = 1:numel(col_idxs)
        labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
    end
    xticklabels(models);
    xtickangle(55);

    t = Neural(row_idx).name;
    if t(end) == 'f'
        t = [t(1:end-2), '_{feedback\_onset}'];
    else
        t = [t(1:end-2), '_{trial\_onset}'];
    end
    ylabel(t);
    if plot_idx == 4
        title('Fisher z-transformed Spearman rank correlation');
    end
end













%% Interesting visualizations
% show select neural / model pairs that I think are interesting, like bar
% plots
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1, 5, 2, 3, 4]; % ROIs to plot at trial onset
neural_idxs = [neural_idxs; neural_idxs + numel(Neural)/2]; % at feedback too
neural_idxs = neural_idxs(:)';

model_idxs = [1:7, 18, 20, 26, 29, 30, 31, 32, 15, 35, 8];

relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
assert(~relative_to_time || ~control_for_time);

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

for i = 1:numel(neural_idxs)
    neural_idx = neural_idxs(i);
    
    % Get fisher z-transformed correlation coefficients
    %
    zs = squeeze(all_subject_zs(neural_idx, model_idxs, :));
    
    if relative_to_time
        time_zs = squeeze(all_subject_zs(neural_idx, 8, :));
        
        % One-sample t-test against the time zs
        %
        [h, ps, ci, stats] = ttest(zs' - time_zs);
        assert(numel(h) == numel(model_idxs));
        
        % plot the z's relative to time
        %
        zs = time_zs - zs;
    else
        % One-sample t-test them against 0
        %
        [h, ps, ci, stats] = ttest(zs');
        assert(numel(h) == numel(model_idxs));
    end
    
    % Compute SEs
    %
    %[sems, means] = wse(zs'); % WSE's
    means = mean(zs, 2);
    sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % ...or just normal SEMs
    assert(numel(sems) == numel(model_idxs));
    
    % Plot
    %
    subplot(numel(neural_idxs) + 1, 1, i);
    
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    set(gca, 'xtick', xs);
    models = {Model(model_idxs).name};
    labels = {};
    for j = 1:numel(model_idxs)
        labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
    end
    xticklabels(labels);
 
    if ~control_for_time
        set(gca, 'ylim', [0 0.405]);
    else
        set(gca, 'ylim', [-0.1 0.03]);
    end
    set(gca, 'ytick', mean(ylim));
    yticklabels({Neural(neural_idx).name});
    %ytickangle(30);
    if i == 1
        title('Fisher z-transformed Spearman rank correlation');
    end
    if i == numel(neural_idxs)
        subplot(numel(neural_idxs) + 1, 1, i + 1);
        h = bar(means * 0, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        set(gca, 'xtick', xs);
        xticklabels(models);
        xtickangle(45);
    end
end












%% Chan et al.-like bar plots for their models
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1:(numel(Neural) / 2)]; % ROIs to plot -- trial_onset
neural_idxs = [1:(numel(Neural) / 2)] + (numel(Neural) / 2); % ROIs to plot -- feedbac_onset
model_idxs = [1:7, 28, 9, 10, 11, 8]; % models to plot -- Chan et al.
relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
assert(~relative_to_time || ~control_for_time);

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

for i = 1:numel(model_idxs)
    model_idx = model_idxs(i);
    
    % Get fisher z-transformed correlation coefficients
    %
    zs = squeeze(all_subject_zs(neural_idxs, model_idx, :));
    
    if relative_to_time
        time_zs = squeeze(all_subject_zs(neural_idxs, 8, :));
        
        % Paired-sample t-test them against the time zs
        %
        [h, ps, ci, stats] = ttest(zs', time_zs');
        assert(numel(h) == numel(neural_idxs));
        
        % plot the z's relative to time
        %
        zs = time_zs - zs;
    else
        % One-sample t-test them against 0
        %
        [h, ps, ci, stats] = ttest(zs');
        assert(numel(h) == numel(neural_idxs));
    end
    
    % Compute SEs
    %
    %[sems, means] = wse(zs'); % WSE's
    means = mean(zs, 2);
    sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % ...or just normal SEMs
    assert(numel(sems) == numel(neural_idxs));
    
    % Plot
    %
    subplot(numel(model_idxs), 1, i);
    
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    set(gca, 'xtick', xs);
    rois = {Neural(neural_idxs).name};
    if i == numel(model_idxs)
        xticklabels(rois);
        xtickangle(45);
    else
        labels = {};
        for j = 1:numel(neural_idxs)
            labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
        end
        xticklabels(labels);
    end
 
    if ~control_for_time
        set(gca, 'ylim', [0 0.405]);
    else
        set(gca, 'ylim', [-0.2 0.05]);
    end
    set(gca, 'ytick', mean(ylim));
    yticklabels({Model(model_idx).name});
    %ytickangle(30);
    if i == 1
        title('Fisher z-transformed Spearman rank correlation');
    end
end



