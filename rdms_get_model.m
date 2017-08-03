function Model = rdms_get_model(data, metadata, which_rows)

% Compute the RDMs for different models
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

% weights before the update (i.e. prior weights) for all structures
%
ww_before = [simulated.ww_before{1} simulated.ww_before{2} simulated.ww_before{3}];
ww_before = ww_before + rand(size(ww_before)) * 0.001; % so cosine metric works
[weightsBeforeRDMs, avgWeightsBeforeRDM] = compute_rdms(ww_before, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsBeforeRDMs;
Model(model_idx).RDM = avgWeightsBeforeRDM;
Model(model_idx).name = 'weightsBefore';
Model(model_idx).color = [0 1 0];



%% models at feedback onset
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
outcomes = repmat(data.outcome, 1, size(simulated.valuess, 2));
PEs = outcomes - simulated.valuess;
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
outcomes = repmat(data.outcome, 1, size(simulated.valuess, 2));
PEs2 = (outcomes - simulated.valuess) .^ 2;
PEs2 = PEs2 + rand(size(PEs2)); % so cosine metric works
[pes2RDMs, avgPes2RDM] = compute_rdms(PEs2, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = pes2RDMs;
Model(model_idx).RDM = avgPes2RDM;
Model(model_idx).name = 'PEsSquared';
Model(model_idx).color = [0 1 0];

% log(PEs^2)
%
outcomes = repmat(data.outcome, 1, size(simulated.valuess, 2));
logPEs2 = log((outcomes - simulated.valuess) .^ 2 + 0.001); % so we don't take log(0)
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

% weights after the update (i.e. posterior weights) for all structures
%
ww_after = [simulated.ww_after{1} simulated.ww_after{2} simulated.ww_after{3}];
ww_after = ww_after + rand(size(ww_after)) * 0.001; % so cosine metric works
[weightsAfterRDMs, avgWeightsAfterRDM] = compute_rdms(ww_after, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsAfterRDMs;
Model(model_idx).RDM = avgWeightsAfterRDM;
Model(model_idx).name = 'weightsAfter';
Model(model_idx).color = [0 1 0];

% prior weights for strucutre corresponding to condition (see GLM 148)
%
max_dim = max([size(simulated.ww_before{1}, 2), size(simulated.ww_before{2}, 2), size(simulated.ww_before{3}, 2)]);
ww_prior = zeros(size(data.condition, 1), max_dim);
for i = 1:size(simulated.ww_before{1}, 1)
    M = current_structure(i);
    ww_prior(i, 1:size(simulated.ww_before{M}, 2)) = simulated.ww_before{M}(i,:);
    % TODO FIXME to make the weights for M1 comparable to M2 and M3, include them twice
    % i.e. when comparing a trial from M1 with a trial from M2, we'll be comparing both weights for x1 and x2 from M1
    % with M2's weights for x1 and x2 in c1 and in c2. I.e. weight vectors
    % are M1 = [x1 x2 x1 x2], M2 = [x1c1 x2c1 x1c2 x2c2]
    % unfortunately for M3, this means we'll be comparing x1 and x2 with c1
    % and c2 because the weight vector for M3 = [x1 x2 c1 c2]
    %
    if M == 1
        assert(size(simulated.ww_before{M}, 2) == 2);
        assert(max_dim == 4);
        ww_prior(i, 3:4) = ww_prior(i, 1:2);
    end
end
ww_prior = ww_prior + rand(size(ww_prior)) * 0.001; % so cosine metric works
[weightsPosteriorRDMs, avgWeightsPosteriorRDM] = compute_rdms(ww_prior, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsPosteriorRDMs;
Model(model_idx).RDM = avgWeightsPosteriorRDM;
Model(model_idx).name = 'weightsPrior';
Model(model_idx).color = [0 1 0];

% posterior weights for strucutre corresponding to condition (see GLM 148)
%
max_dim = max([size(simulated.ww_after{1}, 2), size(simulated.ww_after{2}, 2), size(simulated.ww_after{3}, 2)]);
ww_posterior = zeros(size(data.condition, 1), max_dim);
for i = 1:size(simulated.ww_after{1}, 1)
    M = current_structure(i);
    ww_posterior(i, 1:size(simulated.ww_after{M}, 2)) = simulated.ww_after{M}(i,:);
    if M == 1
        assert(size(simulated.ww_before{M}, 2) == 2);
        assert(max_dim == 4);
        ww_posterior(i, 3:4) = ww_posterior(i, 1:2);
    end
end
ww_posterior = ww_posterior + rand(size(ww_posterior)) * 0.001; % so cosine metric works
[weightsPosteriorRDMs, avgWeightsPosteriorRDM] = compute_rdms(ww_posterior, 'cosine', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = weightsPosteriorRDMs;
Model(model_idx).RDM = avgWeightsPosteriorRDM;
Model(model_idx).name = 'weightsPosterior';
Model(model_idx).color = [0 1 0];


% prior Sigma for strucutre corresponding to condition (see GLM 148)
%
%max_dim = max([numel(simulated.Sigma_before{1}(:,:,1)), numel(simulated.Sigma_before{2}(:,:,1)), numel(simulated.Sigma_before{3}(:,:,1))]);
%Sigma_prior = zeros(size(data.condition, 1), max_dim); % it's ok to pad with zeros when the weight vector is smaller
%for i = 1:size(simulated.Sigma_before{1}, 3)
%    M = current_structure(i);
%    s = simulated.Sigma_before{M}(:,:,i);
%    s = s(:)'; % flatten the covariance matrix
%    Sigma_prior(i, 1:numel(s)) = s;
%end
%Sigma_prior = Sigma_prior + rand(size(Sigma_prior)) * 0.001; % so cosine metric works



disp('Computed model RDMs.');
toc
