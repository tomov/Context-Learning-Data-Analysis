% same as kl_structure_learning.m but for the weights KL divergence
% TODO dedupe NOTE that we flipped subj_idx and run in the kl_betas_roi and test_liks
%

%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Simulate behavior
%
prior_variance = 0.1249;
inv_softmax_temp = 2.0064;
params = [prior_variance inv_softmax_temp];
which_structures = logical([1 1 1 0]);
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows, false);

%% Some random voxels as controls
%
rand_voxels = gen_rand_voxels(fullfile('masks', 'mask.nii'), 20);

%% The KL_weight betas
%

% The peak voxel from each ROI from the contrast (from the results table)
%
KL_weights_rois = {'Temporal_Inf_L', ...
    'Parietal_Inf_L', ...
	'Temporal_Inf_L', ...
	'Supp_Motor_Area_R', ...
	'Frontal_Sup_Medial_L', ...
	'Supp_Motor_Area_L', ...
	'Location not in atlas', ...
	'Location not in atlas', ...
	'Frontal_Inf_Oper_L', ...
	'Frontal_Mid_2_L', ...
	'Frontal_Inf_Tri_R', ...
	'Frontal_Inf_Tri_R', ...
	'Cerebelum_8_R', ...
	'Frontal_Inf_Orb_2_L', ...
	'Pallidum_L', ...
	'Frontal_Inf_Tri_L'};

KL_weights_peak_voxels = [-48    -52    -16; ...
	-32    -56    46; ...
	-42    -34    -28; ...
	10    14    46; ...
	-8    22    42; ...
	-8    10    68; ...
	30    24    0; ...
	16    6    -6; ...
	-42    8    26; ...
	-36    4    58; ...
	36    20    28; ...
	50    34    24; ...
	30    -68    -50; ...
	-46    16    -4; ...
	-20    4    -4; ...
	-44    40    0];

KL_weights_betas = load_run_betas(148, 'KL_weights', [KL_weights_peak_voxels; rand_voxels]);

%% The KL_posterior betas, voxels taken from 'surprise - wrong' contrast
%
KL_posterior_error_rois = {'Angular_R', ...
    'Parietal_Inf_R', ...
    'Frontal_Mid_2_L', ...
    'Location not in atlas', ...
    'Frontal_Mid_2_R', ...
    'OFCmed_R', ...
    'Frontal_Mid_2_R'};

KL_posterior_error_peak_voxels = [34  -68 52; ...
    40  -46 38; ...
    -42 54  2; ...
    -30 62  22; ...
    36  54  0; ...
    18  42  -16; ...
    52  32  32];

KL_posterior_error_betas = load_run_betas(123, 'surprise', [KL_posterior_error_peak_voxels; rand_voxels]);

%% KL_posterior betas for 'surprise' contrast only
%
KL_posterior_rois = { ...
    'Parietal_Sup_R', ...
    'Parietal_Inf_R', ...
    'Location', ...
    'Frontal_Inf_Oper_R', ...
    'Frontal_Mid_2_R', ...
    'Frontal_Mid_2_L', ...
    'Frontal_Sup_2_L', ...
    'Frontal_Inf_Orb_2_L', ...
    'Frontal_Mid_2_R', ...
    'OFCant_R', ...
    'Parietal_Inf_L', ...
    'Parietal_Inf_L', ...
    'Cerebelum_Crus2_L', ...
    'Frontal_Inf_Tri_L', ...
    'Precentral_L', ...
    'Frontal_Mid_2_L', ...
    'Cerebelum_8_L', ...
    'Cerebelum_Crus1_L'};

KL_posterior_peak_voxels = [ ...
    32    -66    50; ...
    42    -46    40; ...
    32    22    0; ...
    48    22    32; ...
    34    12    52; ...
    -40    56    4; ...
    -24    58    -10; ...
    -46    34    -10; ...
    36    54    0; ...
    20    48    -16; ...
    -30    -62    40; ...
    -46    -38    44; ...
    -6    -78    -30; ...
    -38    16    26; ...
    -38    0    40; ...
    -50    34    28; ...
    -32    -68    -54; ...
    -30    -62    -32];

KL_posterior_betas = load_run_betas(123, 'surprise', [KL_posterior_peak_voxels; rand_voxels]);

%% Get the behavioral correlates
%
test_liks = nan(metadata.N, metadata.runsPerSubject);
test_RTs = nan(metadata.N, metadata.runsPerSubject);

% Iterate over subjects
%
subj_idx = 0; % 1..20
for subj = goodSubjects % 1..25
    subject = metadata.allSubjects(subj); % 'con001' ... 'con025'
    subj_trials = data.which_rows & strcmp(data.participant, subject);
    subj_idx = subj_idx + 1;
    fprintf('subj %d (idx %d)\n', subj, subj_idx);

    % Iterate over runs
    %
    for run = 1:metadata.runsPerSubject
        fprintf(' RUN %d\n', run);
        run_trials = subj_trials & data.runId == run;
        condition = data.contextRole(run_trials);
        condition = condition{1};
		assert(sum(run_trials) == metadata.trialsPerRun);

        % Get the test trial likelihood
        %
        run_test_trials = run_trials & ~data.isTrain & ~data.timeout;
        
        X_fixed = data.chose_sick(run_test_trials); % actual subject choice on each trial
        P = simulated.pred(run_test_trials); % probability of subject choice on each trial
        assert(numel(X_fixed) <= metadata.testTrialsPerRun);
        
        liks = binopdf(X_fixed, 1, P);
        assert(numel(liks) == numel(X_fixed));
        avg_loglik = mean(log(liks)); % average to account for timeouts
        
        test_liks(subj_idx, run) = avg_loglik; 
        
        % Get the test trial RTs
        %
        test_RTs(subj_idx, run) = mean(data.response.rt(run_test_trials));
    end
end


save('results/kl_weights.mat');

%% within-subject analysis using a linear mixed effects model and/or t-tests
%

load('results/kl_weights.mat');

% pick which one to analyze TODO refactor
%
analyze_KL_weights_and_not_KL_posterior = false; % #KNOB

if analyze_KL_weights_and_not_KL_posterior
    KL_betas = KL_weights_betas;
    rois = KL_weights_rois;
else
    KL_betas = KL_posterior_betas;
    rois = KL_posterior_rois;
end

% Do the correlations
%

lmes = {};
lme_means = [];
lme_sems = [];
lme_ps = [];

glm_means = [];
glm_sems = [];
glm_ps = [];

ttest_means = [];
ttest_sems = [];
ttest_ps = [];

for roi = 1:size(KL_betas, 3)
    kl_betas_roi = KL_betas(:, :, roi);

    lme_liks = [];
    lme_betas = [];
    lme_ids = [];
    
    fisher_rs = [];

    for subj_idx = 1:metadata.N
        x = test_liks(subj_idx, :)';
        y = kl_betas_roi(subj_idx, :)';
        
        %x = zscore(x);
        %y = zscore(y);

        % LME model within-subject analysis
        %
        lme_liks = [lme_liks; x];
        lme_betas = [lme_betas; y];
        lme_ids = [lme_ids; ones(size(x,1),1) * subj_idx];
        
        % Simple within-subject t-test
        %
        r = corrcoef(x, y);
        r = r(1,2);
        fisher_r = atanh(r);
        fisher_rs = [fisher_rs; fisher_r];
    end

    % LME model
    %
    tbl = array2table([lme_betas, lme_liks, lme_ids], 'VariableNames', {'Beta', 'Likelihood', 'Subject'});
    %formula = 'Beta ~ Likelihood + (Likelihood|Subject)';
    formula = 'Beta ~ Likelihood + (Likelihood|Subject)';
    lme = fitlme(tbl, formula);
    lmes{roi} = lme;

    p = lme.Coefficients(2, 'pValue').pValue;
    coef = lme.Coefficients(2, 'Estimate').Estimate;
    lme_means = [lme_means; coef];
    lme_sems = [lme_sems; (lme.Coefficients(2, 'Upper').Upper - lme.Coefficients(2, 'Lower').Lower) / 2];
    lme_ps = [lme_ps; p];
    
    % t-test on fisher z-transformed correlation coefficients
    %
    [h, p, ci, stats] = ttest(fisher_rs);
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(fisher_rs)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
    
    % manually run LME with the GLM function
    % WRONG -- betas for random effects come from normal distribution
    %
    %{
    o = ones(size(lme_betas, 1), 1);
    X_fixed = [o lme_liks];
    X = X_fixed;
    for subj_idx = 1:metadata.N
        X_subj = X_fixed .* (lme_ids == subj_idx);
        X = [X X_subj];
    end
    [b, dev, stats] = glmfit(X, lme_betas, 'normal', 'constant', 'off');
    glm_means = [glm_means; b(2)];
    glm_sems = [glm_sems; stats.se(2)];
    glm_ps = [glm_ps; stats.p(2)];
    %}
end

%
% LME plots
%

figure;

% plot correlation coefficients for LME
%
subplot(2, 1, 1);

h = bar(lme_means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
errorbar(xs, lme_means, lme_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;

set(gca, 'xtick', xs);
ylabel('Fixed effect');
xticklabels([strrep(rois, '_', '\_'), repmat({'random'}, 1, size(rand_voxels, 1))]);
xtickangle(60);
xlabel('voxel');
title('KL_weight betas correlated with test log likelihood: LME analysis', 'Interpreter', 'none');

% plot p-values
%
subplot(2, 1, 2);

h = bar(lme_ps, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
plot([0 max(xs) + 1],[0.05 0.05],'k--');
hold off;

set(gca, 'xtick', xs);
ylabel('p-value');
xticklabels([strrep(rois, '_', '\_'), repmat({'random'}, 1, size(rand_voxels, 1))]);
xtickangle(60);
xlabel('voxel');

% plot random effects for first ROI
%
%{
subplot(3, 1, 3);

roi_idx = 1;
lme = lmes{roi_idx};

[randMeans, ~, randStats] = randomEffects(lme);
assert(numel(randMeans) == metadata.N * 2);
randSems = (randStats(:, 'Upper').Upper - randStats(:, 'Lower').Lower) / 2;

h = bar(reshape(randMeans, [2 20])', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
xs = sort([xs - 0.2, xs + 0.2])';
hold on;
errorbar(xs, randMeans, randSems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;

set(gca, 'xtick', h(1).XData);
ylabel('Random effect');
xtickangle(60);
xlabel('subject');

[fixedMeans, ~, fixedStats] = fixedEffects(lme);
%}


%
% t-test plots
%

figure;

% plot correlation coefficients for t-tests
%
subplot(2, 1, 1);

h = bar(ttest_means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
errorbar(xs, ttest_means, ttest_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;

set(gca, 'xtick', xs);
ylabel('Fixed effect');
xticklabels([strrep(rois, '_', '\_'), repmat({'random'}, 1, size(rand_voxels, 1))]);
xtickangle(60);
xlabel('voxel');
title('KL_weight betas correlated with test log likelihood: t-test', 'Interpreter', 'none');

% plot p-values
%
subplot(2, 1, 2);

h = bar(ttest_ps, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
plot([0 max(xs) + 1],[0.05 0.05],'k--');
hold off;

set(gca, 'xtick', xs);
ylabel('p-value');
xticklabels([strrep(rois, '_', '\_'), repmat({'random'}, 1, size(rand_voxels, 1))]);
xtickangle(60);
xlabel('voxel');



