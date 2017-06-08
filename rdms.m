% Compare RDMs from different ROIs with model RDMs
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
Neural(neural_idx).RDMs = hippocampusRDMs;
Neural(neural_idx).RDM = avgHippocampusRDM;
Neural(neural_idx).name = 'hippocampus';
Neural(neural_idx).color = [0 1 0];

ofc_mask = load_mask('masks/ofc.nii');
ofc_betas = get_betas_submask(ofc_mask, whole_brain_trial_onset_betas);
[ofcRDMs, avgOfcRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = ofcRDMs;
Neural(neural_idx).RDM = avgOfcRDM;
Neural(neural_idx).name = 'OFC';
Neural(neural_idx).color = [0 1 0];

med_ofc_mask = load_mask('masks/med_ofc.nii');
med_ofc_betas = get_betas_submask(med_ofc_mask, whole_brain_trial_onset_betas);
[medOfcRDMs, avgMedOfcRDM] = compute_rdms(med_ofc_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = medOfcRDMs;
Neural(neural_idx).RDM = avgMedOfcRDM;
Neural(neural_idx).name = 'mOFC';
Neural(neural_idx).color = [0 1 0];

rectus_mask = load_mask('masks/rectus.nii');
rectus_betas = get_betas_submask(rectus_mask, whole_brain_trial_onset_betas);
[rectusRDMs, avgRectusRDM] = compute_rdms(rectus_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = rectusRDMs;
Neural(neural_idx).RDM = avgRectusRDM;
Neural(neural_idx).name = 'Rectus';
Neural(neural_idx).color = [0 1 0];

vmpfc_mask = load_mask('masks/vmpfc.nii');
vmpfc_betas = get_betas_submask(vmpfc_mask, whole_brain_trial_onset_betas);
[vmpfcRDMs, avgVmpfcRDM] = compute_rdms(vmpfc_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = vmpfcRDMs;
Neural(neural_idx).RDM = avgVmpfcRDM;
Neural(neural_idx).name = 'vmPFC';
Neural(neural_idx).color = [0 1 0];

striatum_mask = load_mask('masks/striatum.nii');
striatum_betas = get_betas_submask(striatum_mask, whole_brain_trial_onset_betas);
[striatumRDMs, avgStriatumRDM] = compute_rdms(striatum_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = striatumRDMs;
Neural(neural_idx).RDM = avgStriatumRDM;
Neural(neural_idx).name = 'Striatum';
Neural(neural_idx).color = [0 1 0];

pallidum_mask = load_mask('masks/pallidum.nii');
pallidum_betas = get_betas_submask(pallidum_mask, whole_brain_trial_onset_betas);
[pallidumRDMs, avgPallidumRDM] = compute_rdms(pallidum_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = pallidumRDMs;
Neural(neural_idx).RDM = avgPallidumRDM;
Neural(neural_idx).name = 'Pallidum';
Neural(neural_idx).color = [0 1 0];

v1_mask = load_mask('masks/v1.nii');
v1_betas = get_betas_submask(v1_mask, whole_brain_trial_onset_betas);
[v1RDMs, avgV1RDM] = compute_rdms(v1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = v1RDMs;
Neural(neural_idx).RDM = avgV1RDM;
Neural(neural_idx).name = 'V1';
Neural(neural_idx).color = [0 1 0];

m1_mask = load_mask('masks/m1.nii');
m1_betas = get_betas_submask(m1_mask, whole_brain_trial_onset_betas);
[m1RDMs, avgM1RDM] = compute_rdms(m1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = m1RDMs;
Neural(neural_idx).RDM = avgM1RDM;
Neural(neural_idx).name = 'M1';
Neural(neural_idx).color = [0 1 0];

s1_mask = load_mask('masks/s1.nii');
s1_betas = get_betas_submask(s1_mask, whole_brain_trial_onset_betas);
[s1RDMs, avgS1RDM] = compute_rdms(s1_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = s1RDMs;
Neural(neural_idx).RDM = avgS1RDM;
Neural(neural_idx).name = 'S1';
Neural(neural_idx).color = [0 1 0];

fusiform_mask = load_mask('masks/fusiform.nii');
fusiform_betas = get_betas_submask(fusiform_mask, whole_brain_trial_onset_betas);
[fusiformRDMs, avgFusiformRDM] = compute_rdms(fusiform_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = fusiformRDMs;
Neural(neural_idx).RDM = avgFusiformRDM;
Neural(neural_idx).name = 'Fusiform';
Neural(neural_idx).color = [0 1 0];

angular_mask = load_mask('masks/angular.nii');
angular_betas = get_betas_submask(angular_mask, whole_brain_trial_onset_betas);
[angularRDMs, avgAngularRDM] = compute_rdms(angular_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = angularRDMs;
Neural(neural_idx).RDM = avgAngularRDM;
Neural(neural_idx).name = 'Angular';
Neural(neural_idx).color = [0 1 0];

mid_front_mask = load_mask('masks/mid_front.nii');
mid_front_betas = get_betas_submask(mid_front_mask, whole_brain_trial_onset_betas);
[midFrontRDMs, avgMidFrontRDM] = compute_rdms(mid_front_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = midFrontRDMs;
Neural(neural_idx).RDM = avgMidFrontRDM;
Neural(neural_idx).name = 'MidFront';
Neural(neural_idx).color = [0 1 0];

dl_sup_front_mask = load_mask('masks/dl_sup_front.nii');
dl_sup_front_betas = get_betas_submask(dl_sup_front_mask, whole_brain_trial_onset_betas);
[dlSupFrontRDMs, avgDlSupFrontRDM] = compute_rdms(dl_sup_front_betas, 'cosine', data, metadata, which_rows);
neural_idx = neural_idx + 1;
Neural(neural_idx).RDMs = dlSupFrontRDMs;
Neural(neural_idx).RDM = avgDlSupFrontRDM;
Neural(neural_idx).name = 'dlSupFront';
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

% D_KL: -abs(D1 - D2) 
%
[surpriseRDMs, avgSurpriseRDM] = compute_rdms(simulated.surprise, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = surpriseRDMs;
Model(model_idx).RDM = avgSurpriseRDM;
Model(model_idx).name = 'KL';
Model(model_idx).color = [0 1 0];

% log of D_KL: -abs(D1 - D2)
%
logSurprise = log(simulated.surprise + 0.001); % so we don't take log(0)
[logSurpriseRDMs, avgLogSurpriseRDM] = compute_rdms(logSurprise, 'euclidean', data, metadata, which_rows);
model_idx = model_idx + 1;
Model(model_idx).RDMs = logSurpriseRDMs;
Model(model_idx).RDM = avgLogSurpriseRDM;
Model(model_idx).name = 'logKL';
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


%% Within-subject comparison
% compare each neural and model RDMs for subject separately using
% Spearman's rank coefficient, then find which ones are significant
%
trig = logical(triu(ones(20,20), 1)); % upper right triangle, excluding diagonal

table_Rho = []; % average Spearman's rho for each ROI for each model
table_H = []; % result of hypothesis test for each ROI for each model -- is the correlation significant?
table_P = []; % p-value of hypothesis test for each ROI for each model
for neural_idx = 1:numel(Neural)

    all_rhos = []; % Spearman rhos: row = model, col = subject
    for model_idx = 1:numel(Model)

        % Compute a Spearman's rank correlation for each subject separately
        %
        rhos = [];
        for subj = 1:metadata.N
            model_RDM = Model(model_idx).RDMs(:,:,subj);
            neural_RDM = Neural(neural_idx).RDMs(:,:,subj);

            x = model_RDM(trig);
            y = neural_RDM(trig);
            rho = corr(x, y, 'type', 'Spearman');
            rhos = [rhos, rho];

            % or for each run even?? how to you get WSE's otherwise?
            %{
            for run = 1:metadata.runsPerSubject
                s = (run - 1) * metadata.trainingTrialsPerRun + 1;
                e = run * metadata.trainingTrialsPerRun;
                model_subRDM = Model(model_idx).RDMs(s:e, s:e, subj);
                neural_subRDM = Neural(neural_idx).RDMs(s:e, s:e, subj);
            end
            %}
        end
        all_rhos = [all_rhos; rhos];
    end
    
    % Group-level analysis
    %
    fisher_all_rhos = atanh(all_rhos);
    [h, ps, ci, stats] = ttest(fisher_all_rhos');

    table_Rho = [table_Rho; mean(fisher_all_rhos')];
    table_H = [table_H; h];
    table_P = [table_P; ps];
end


Rho = array2table(table_Rho, 'RowNames', {Neural.name}, 'VariableNames', {Model.name});
H = array2table(table_H, 'RowNames', {Neural.name}, 'VariableNames', {Model.name});
P = array2table(table_P, 'RowNames', {Neural.name}, 'VariableNames', {Model.name});
