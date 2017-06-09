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

% values -- normalized correlation (cosine) is not a good measure here
%
[valsRDMs, avgValsRDM] = compute_rdms(simulated.valuess(:, which_structures), 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = valsRDMs;
Model(model_idx).RDM = avgValsRDM;
Model(model_idx).name = 'values';
Model(model_idx).color = [0 1 0];

% weights -- normalized correlation (cosine) is not a good measure here
%
ww = [simulated.ww1 simulated.ww2 simulated.ww3];
[weightsRDMs, avgWeightsRDM] = compute_rdms(ww, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsRDMs;
Model(model_idx).RDM = avgWeightsRDM;
Model(model_idx).name = 'weights';
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

% PEs -- normalized correlation (cosine) is not a good measure here
%
PEs = data.outcome - simulated.valuess;
[pesRDMs, avgPesRDM] = compute_rdms(PEs, 'euclidean', data, metadata, which_rows);
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

% PEs^2 -- normalized correlation (cosine) is not a good measure here
%
PEs2 = (data.outcome - simulated.valuess) .^ 2;
[pes2RDMs, avgPes2RDM] = compute_rdms(PEs2, 'euclidean', data, metadata, which_rows);
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
Model(model_idx).name = 'new_value';
Model(model_idx).color = [0 1 0];

% New values -- normalized correlation (cosine) is not a good measure here
%
[newValsRDMs, avgNewValsRDM] = compute_rdms(simulated.new_valuess(:, which_structures), 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = newValsRDMs;
Model(model_idx).RDM = avgNewValsRDM;
Model(model_idx).name = 'new_values';
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



% show the model RDMs
% 
showRDMs(Model, 2);

%
% Second-order similarity matrix from the RDMs
% compare neural RDMs with the model RDMs
%

%% Compare the average RDMs only
%
userOptions.RDMcorrelationType= 'Spearman';
userOptions.analysisName = 'blah';
userOptions.rootPath = '~/Downloads/'; % TODO how to turn off saving the figure?
corrMat = pairwiseCorrelateRDMs({Neural, Model}, userOptions, struct('figureNumber', 3,'fileName',[]));

%% Within-subject / Between-subject RDM comparison
% compare each neural and model RDMs for subject separately using
% Spearman's rank coefficient, then find which ones are significant
%
tic;

within_subject = true; % compute correlation for each run separately i.e. no comparisons of representations across runs
all_vs_all = false; % compare neural vs. neural and model vs. model representations too

trig = logical(triu(ones(metadata.runsPerSubject * metadata.trainingTrialsPerRun), 1)); % upper right triangle, excluding diagonal, for all runs
run_trig = logical(triu(ones(metadata.trainingTrialsPerRun), 1)); % upper right triangle, excluding diagonal, for a single run

table_Rho = []; % average Spearman's rho for each ROI for each model
table_H = []; % result of hypothesis test for each ROI for each model -- is the correlation significant?
table_P = []; % p-value of hypothesis test for each ROI for each model

between_subject_rhos = nan(numel(Neural), numel(Model), metadata.N);
within_subject_rhos = nan(numel(Neural), numel(Model), metadata.N, metadata.runsPerSubject);

if all_vs_all
    all = [Neural, Model];
    rows = all;
    cols = all;
else
    rows = Neural;
    cols = Model;
end

for row_idx = 1:numel(rows)

    models_subjs_rhos = []; % Spearman rhos: row = model, col = subject
    for col_idx = 1:numel(cols)

        % Compute a Spearman's rank correlation for each subject separately
        %
        subjs_rhos = []; % Spearman rhos: col = subject, one rho per
        for subj = 1:metadata.N
            row_RDM = rows(row_idx).RDMs(:,:,subj);
            col_RDM = cols(col_idx).RDMs(:,:,subj);
            assert(isequal(size(row_RDM), size(col_RDM)));
            assert(isequal(size(trig), size(row_RDM)));

            if ~within_subject
                % Look at RDMs for all runs simultaneously
                %
                x = row_RDM(trig);
                y = col_RDM(trig);
                subj_rho = corr(x, y, 'type', 'Spearman');
                subjs_rhos = [subjs_rhos, subj_rho];
                
                between_subject_rhos(row_idx, col_idx, subj) = subj_rho;
            else
                % Look at RDM for each run separately
                %
                % or for each run even?? how to you get WSE's otherwise?
                % but it doesn't make sense for us b/c we have 1 condition per
                % run => the RDMs will look similar across runs, even if the
                % posteriors are actually different
                %
                runs_rhos = []; % one rho per run
                for run = 1:metadata.runsPerSubject
                    s = (run - 1) * metadata.trainingTrialsPerRun + 1;
                    e = run * metadata.trainingTrialsPerRun;
                    row_subRDM = rows(row_idx).RDMs(s:e, s:e, subj);
                    col_subRDM = cols(col_idx).RDMs(s:e, s:e, subj);
                    assert(isequal(size(row_subRDM), size(col_subRDM)));
                    assert(isequal(size(run_trig), size(row_subRDM)));

                    x = row_subRDM(run_trig);
                    y = col_subRDM(run_trig);
                    run_rho = corr(x, y, 'type', 'Spearman');
                    runs_rhos = [runs_rhos, run_rho];
                    
                    within_subject_rhos(row_idx, col_idx, subj, run) = run_rho;
                end
                subj_rho = mean(runs_rhos); % subject rho = mean rho across runs TODO use WSE
                subjs_rhos = [subjs_rhos, subj_rho];
            end
        end
        models_subjs_rhos = [models_subjs_rhos; subjs_rhos];
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

%% Visualize analysis
%

tabs = {table_Rho, table_P};
titles = {'Representational similarity match (Spearman rank correlation)', 'P-value (one-sample t-test)'};

for i = 1:numel(tabs)
    figure;
    if i == 2
        t = tabs{i};
        if all_vs_all
            imagesc(log10(t), [-18 0.1]);
        else
            imagesc(log10(t));
        end    
        c = colorbar;
        y = get(c, 'Ytick');
        set(c, 'YTickLabel', arrayfun(@(x) ['10\^', num2str(x)], y, 'UniformOutput', false));
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

    yticklabels({rows.name})
    set(gca, 'ytick', 1:numel({rows.name}));
    ylabel('ROI_{event}, t = trial onset, f = feedback onset');
    
    title(titles{i});
end


%% More visualization for within-subject analysis only
%
%{
assert(within_subject);

neural_idx = 18; % striatum_f
w = atanh(within_subject_rhos);
w = w - w(neural_idx,:,:,:);

model_idx = 15;
for ni = 1:numel(Neural)
    w(ni, model_idx, :, :)
end
%}