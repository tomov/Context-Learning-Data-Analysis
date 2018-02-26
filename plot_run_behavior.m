function plot_run_behavior(subj, run)
% Plot behavioral data for a single run; hybrid of plot_spm_regressors and plot_behavior
% for sanity checks; run before & compare with plot_spm_regressors (before running full GLM)
%

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
which_rows = data.which_rows & strcmp(metadata.allSubjects{subj}, data.participant) & data.runId == run;


figure;

subplot_idx = 0;


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

% show stimulus sequence
% this is ConLearn-specific 
% TODO dedupe with plot_spm_regressors
%

hold on;

which_train = ~data.drop & data.isTrain & strcmp(data.participant, metadata.allSubjects{subj}) & data.runId == run;
condition = data.contextRole(which_train);
condition = condition{1};

is_sick = strcmp(data.sick(which_rows), 'Yes');
is_train = which_train(which_rows);
resps = data.response.keys(which_rows);
chose_sick = double(strcmp(resps, 'left'));
is_timeout = strcmp(resps, 'None');
corr = data.response.corr(which_rows);

onsets = cellfun(@str2num, data.actualChoiceOnset(which_rows)');
for t=onsets
    plot([t t], [-2 2], '--', 'Color', [0.8 0.8 0.8]);
end

feedback_onsets = cellfun(@str2num, data.actualFeedbackOnset(which_train)');
for i=1:length(feedback_onsets)
	if corr(i)
        color = 'blue';
    else
        color = 'red';
    end
    plot([feedback_onsets(i) feedback_onsets(i)], [-2 2], 'Color', color);
end

stims = strcat('x', num2str(data.cueId(which_rows) + 1), 'c', num2str(data.contextId(which_rows) + 1));
show_text_nontimeout = @(which_trials, color) text(onsets(which_trials & ~is_timeout), chose_sick(which_trials & ~is_timeout), stims(which_trials & ~is_timeout,:), 'Color', color);
show_text_timeout = @(which_trials, color) text(onsets(which_trials & is_timeout), 0.5 * ones(sum(which_trials & is_timeout), 1), stims(which_trials & is_timeout,:), 'Color', color);

% sick training trials
show_text_nontimeout(is_sick & is_train, [0 0.5 0]);
show_text_timeout(is_sick & is_train, [0 0.5 0]);
% not sick training trials
show_text_nontimeout(~is_sick & is_train, [0.5 0.5 0]);
show_text_timeout(~is_sick & is_train, [0.5 0.5 0]);
% test trials
show_text_nontimeout(~is_train, [0.3 0.3 0.3]);
show_text_timeout(~is_train, [0.3 0.3 0.3]);

ylim([-0.3 1.1]);
yticks([0 0.5 1]);
yticklabels({'chose not sick', 'timeout', 'chose sick'});

% for legend
h = [
     plot(NaN,NaN, '--', 'Color', [0.8 0.8 0.8],'LineWidth', 1); ...
     plot(NaN,NaN, 'Color', 'red','LineWidth', 1); ...
     plot(NaN,NaN, 'Color', 'blue','LineWidth', 1); ...
     plot(NaN,NaN, 'Color', [0 0.5 0],'LineWidth', 2); ...
     plot(NaN,NaN,'Color', [0.5 0.5 0],'LineWidth', 2); ...
     plot(NaN,NaN,'Color', [0.3 0.3 0.3],'LineWidth', 2)];
legend(h, {'trial onset', 'WRONG/TIMEOUT feedback', 'CORRECT feedback', 'sick outcome', 'not sick outcome', 'test trial'});

title(sprintf('subject %d, run %d, condition: %s', subj, run, upper(condition)));
hold off;



%
% GLM 154
%

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

[data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results.mat'), 1, [1 1 1 0 ], which_rows);

P = simulated.P(which_train,:);
plot(feedback_onsets, P, 'o-', 'LineWidth', 2);
ylabel('P_n');
title('M1, M2, M3: Posterior');
legend({'M1', 'M2', 'M3', 'M1''', 'M2'''});

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL = simulated.surprise(which_train);
plot(feedback_onsets, KL, 'o-', 'LineWidth', 2);
ylabel('KL(P_n || P_{n-1})');
title('KL structures');


%
% GLM 156
%

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

[data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_reviewer2.mat'), 1, [1 1 0 1 0], which_rows);

P = simulated.P(which_train,:);
plot(feedback_onsets, P, 'o-', 'LineWidth', 2);
ylabel('P_n');
title('M1, M2, M1'': Posterior');
legend({'M1', 'M2', 'M3', 'M1''', 'M2'''});

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL = simulated.surprise(which_train);
plot(feedback_onsets, KL, 'o-', 'LineWidth', 2);
ylabel('KL(P_n || P_{n-1})');
title('KL structures');


%
% GLM 157 & 158
%

[data, metadata, simulated] = simulate_subjects_helper(true, 'results/fit_params_results_simple_collins_5nstarts.mat', 1, 'simple_collins');

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc_given_c = simulated.surprise_Zc_given_c(which_train);
KL_Zs_given_s = simulated.surprise_Zs_given_s(which_train);
plot(feedback_onsets, [KL_Zc_given_c KL_Zs_given_s], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for cluster conditionals P(Z_c|c)');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

PE = simulated.PEs(which_train);
plot(feedback_onsets, PE, 'o-', 'LineWidth', 2);
ylabel('PE');
title('prediction error');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc_and_C = simulated.surprise_Zc_and_C(which_train);
KL_Zs_and_S = simulated.surprise_Zs_and_S(which_train);
plot(feedback_onsets, [KL_Zc_and_C KL_Zs_and_S], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for joint P(Z_c,C)');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc = simulated.surprise_Zc(which_train);
KL_Zs = simulated.surprise_Zs(which_train);
plot(feedback_onsets, [KL_Zc KL_Zs], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for cluster marginal P(Z_c)');

%
% GLM 159
%

[data, metadata, simulated] = simulate_subjects_helper(true, 'results/fit_params_results_simple_collins_25nstarts.mat', 1, 'simple_collins');

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc_given_c = simulated.surprise_Zc_given_c(which_train);
KL_Zs_given_s = simulated.surprise_Zs_given_s(which_train);
plot(feedback_onsets, [KL_Zc_given_c KL_Zs_given_s], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for cluster conditionals P(Z_c|c)');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

PE = simulated.PEs(which_train);
plot(feedback_onsets, PE, 'o-', 'LineWidth', 2);
ylabel('PE');
title('prediction error');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc_and_C = simulated.surprise_Zc_and_C(which_train);
KL_Zs_and_S = simulated.surprise_Zs_and_S(which_train);
plot(feedback_onsets, [KL_Zc_and_C KL_Zs_and_S], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for joint P(Z_c,C)');


subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

KL_Zc = simulated.surprise_Zc(which_train);
KL_Zs = simulated.surprise_Zs(which_train);
plot(feedback_onsets, [KL_Zc KL_Zs], 'o-', 'LineWidth', 2);
legend({'contexts', 'cues'});
title('Collins: KL for cluster marginal P(Z_c)');

%
% GLM 160
%

[data, metadata, simulated] = simulate_subjects_helper(true, 'results/fit_params_results_flat_collins_5nstarts.mat', 1, 'flat_collins');

subplot_idx = subplot_idx + 1;
subplot(14,1,subplot_idx);

PEs = simulated.PEs(which_train);
plot(feedback_onsets, PEs, 'o-', 'LineWidth', 2);
title('Flat: PEs');

