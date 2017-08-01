function [ttest_means, ttest_sems, ttest_ps, ttest_ts, lmes, lme_means, lme_sems, lme_ps] = correlate_neural_and_behavior(neural_activations, neural_names, behavioral_measure, plot_title)
% Within-subject correlation of some measure of neural activity (e.g. KL
% betas) and some measure of behavior (e.g. test choice log likelihood).
% First correlates them for each subject (9 pairs of measures per subject, one per run)
% Then t-tests the Fisher z-transformed correlation coefficients against 0
% (20 z's, one per subject).
%
% INPUT:
% neural_activations = N x R x V vector of neural correlates, where N = # subjects, R = # runs
%                      per subject, and V = # of ROIs
%                      e.g. the output of load_run_betas()
% neural_names = cell array of size V with the name of each ROI
% behavioral_measure = N x R matrix of some behavioral measure
%                      e.g. the output of get_test_behavior()
% plot_title = title of the plot; if empty, no plot is plotted
%
% OUTPUT:
% ttest_means, ttest_sems, ttest_ps, ttest_ts, lmes, lme_means, lme_sems,
% lme_ps = self-explanatory

assert(size(neural_activations, 1) == size(behavioral_measure, 1));
assert(size(neural_activations, 2) == size(behavioral_measure, 2));
assert(size(neural_activations, 3) == numel(neural_names));

% Load behavior
%
[~, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();


% Do the correlations
%

lmes = {};
lme_means = [];
lme_sems = [];
lme_ps = [];

ttest_means = [];
ttest_sems = [];
ttest_ps = [];
ttest_ts = [];

all_fisher_rs = {};

for roi = 1:size(neural_activations, 3)
    kl_betas_roi = neural_activations(:, :, roi);

    lme_liks = [];
    lme_betas = [];
    lme_ids = [];
    
    fisher_rs = [];

    for subj_idx = 1:metadata.N
        x = behavioral_measure(subj_idx, :)';
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
    ttest_ts = [ttest_ts; stats.tstat];
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(fisher_rs)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
    all_fisher_rs{roi} = fisher_rs;
end

%
% LME plots
%
%{
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
title([plot_what, ' betas correlated with test log likelihood: LME analysis'], 'Interpreter', 'none');

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
%}

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

if plot_title
    
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
    xticklabels(strrep(neural_names, '_', '\_'));
    xtickangle(60);
    xlabel('voxel');
    title(plot_title, 'Interpreter', 'none');

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
    xticklabels(strrep(neural_names, '_', '\_'));
    xtickangle(60);
    xlabel('voxel');

end



