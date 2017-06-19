% same as kl_structure_learning.m but for the weights KL divergence
% TODO dedupe NOTE that we flipped subj_idx and run in the kl_betas_roi and test_liks
%

% GLM 148 = condition-specific KL_weights & wrong
%
KL_glm = 148;

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

%% sanity check max voxel from 'KL_weights' contrast
% to make sure our method of extracting these is correct
%
EXPT = context_expt();
modeldir = EXPT.modeldir;
V = spm_vol(fullfile(modeldir, ['model', num2str(KL_glm)], ['con6'], 'spmT_0001.nii')); % T-value map
Y = spm_read_vols(V);
cor = mni2cor([-48    -52    -16], V.mat);
Y(cor(1), cor(2), cor(3)); % sanity check -- should be as seen in ccnl_view Show Results Table
assert(abs(Y(cor(1), cor(2), cor(3)) - 9.032) < 1e-3);

%% Pick the voxels to look at
%

% Some random voxels as controls
%
rand_voxels = gen_rand_voxels(fullfile('masks', 'mask.nii'), 20);

% The peak voxel from each ROI from the contrast (from the results table)
%
rois = {'Temporal_Inf_L', ...
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

peak_voxels = [-48    -52    -16; ...
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

%% Load the KL betas (takes a while)
%
kl_betas = load_run_betas(148, 'KL_weights', [peak_voxels; rand_voxels]);

%% !!!!! use the posterior KL !!!!! TODO refactor
%rois = {'Angular_R', 'Parietal_Inf_R', 'Frontal_Mid_2_L', 'Location not in atlas', 'Frontal_Mid_2_R', 'OFCmed_R', 'Frontal_Mid_2_R'};
%peak_voxels = [34  -68 52; 40  -46 38; -42 54  2; -30 62  22; 36  54  0; 18  42  -16; 52  32  32];
%kl_betas = load_run_betas(123, 'surprise', [peak_voxels; rand_voxels]);

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
        
        X = data.chose_sick(run_test_trials); % actual subject choice on each trial
        P = simulated.pred(run_test_trials); % probability of subject choice on each trial
        assert(numel(X) <= metadata.testTrialsPerRun);
        
        liks = binopdf(X, 1, P);
        assert(numel(liks) == numel(X));
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

lmes = {};
means = [];
sems = [];
ps = [];

ttest_ps = [];
ttest_means = [];
ttest_sems = [];

for roi = 1:size(kl_betas, 3)
    kl_betas_roi = kl_betas(:, :, roi);

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
    formula = 'Beta ~ Likelihood + (Likelihood|Subject)';
    lme = fitlme(tbl, formula);
    lmes{roi} = lme;

    p = lme.Coefficients(2, 'pValue').pValue;
    coef = lme.Coefficients(2, 'Estimate').Estimate;
    means = [means; coef];
    sems = [sems; (lme.Coefficients(2, 'Upper').Upper - lme.Coefficients(2, 'Lower').Lower) / 2];
    ps = [ps; p];
    
    % t-test on fisher z-transformed correlation coefficients
    %
    [h, p, ci, stats] = ttest(fisher_rs);
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(fisher_rs)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
end

%
% LME plots
%

figure;

% plot correlation coefficients for LME
%
subplot(2, 1, 1);

h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
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

h = bar(ps, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
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
