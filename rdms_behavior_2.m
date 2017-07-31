% TODO dedupe with rdms_behavior.m and rdms_second_order.m
%
% Idea: RDMs of ROIs vs. models (e.g. posterior) on a run-by-run basis;
% correlate with behavior 
% same as kl_structure_learning effect but multivoxel pattern
%

close all;
clear all;

utils;

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

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


predict = @(V_n) softmax(V_n, params(2));


simulated = simulate_subjects(data, metadata, params, which_structures);


%% Get the neural RDMs
%
clear masks;
mask_idx = 0;

% GLM 123 (KL_structures) clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_1172voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-cluster-1172';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_458voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-cluster-458';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_192voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-cluster-192';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_59voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-cluster-49';

% searchlight -- posterior @ feedback_onset clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_911voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-cluster-911';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_290voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-cluster-290';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_60voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-cluster-60';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_1989voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-dlPFC-cluster-1989';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_117voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-dlPFC-cluster-117';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-50_y=24_z=26_1038voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-dlPFC-cluster-1038';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-52_y=22_z=26_61voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-dlPFC-cluster-61';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=4_y=22_z=44_253voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-ACC-cluster-253';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=6_y=22_z=44_37voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-ACC-cluster-37';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=36_602voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-cluster-602';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_119voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-cluster-119';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_10voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-cluster-10';


% searchlight -- prior @ trial_onset clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-54_y=12_z=26_315voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-cluster-315';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=14_z=22_28voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-cluster-28';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=12_z=22_279voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-cluster-279';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC1_ClusterMask_searchlight_tmap_x=-36_y=10_z=24_17voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC1-cluster-17';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC2_ClusterMask_searchlight_tmap_x=-56_y=14_z=26_34voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC2-cluster-34';


%Neural = rdms_get_neural_2(masks, data, metadata, which_trials, false, false);
Neural = rdms_get_neural(data, metadata, which_trials, false, false); % use t-maps, use nosmooth
%Neural = [Neural, Neural_controls];
showRDMs(Neural, 1);

%% Get the model RDMs
%
Model = rdms_get_model(data, metadata, which_trials);
showRDMs(Model, 2);

control_model_idxs = [8, 12]; % #KNOB control for time and run
assert(isequal(Model(8).name, 'time'));
assert(isequal(Model(12).name, 'run'));

model_idx = 1;
assert(isequal(Model(1).name, 'posterior'));


%% find how similar representations are to posterior (model_idx = 1) in different ROIs per-run
%
rhos = nan(metadata.runsPerSubject, metadata.N, numel(Neural));

goodSubjects = getGoodSubjects();
subjs = metadata.allSubjects(goodSubjects);

% t1, t2 are trial indices corresponding to the pairs of trials in each
% cell of the RDM (for a single subject)
% used to quickly generate RDM masks and extract sub-RDMs.
% Assumes the same trials were used from all subjects (so it takes subj 1
% for convenience)
%
which_trials_per_subj = which_trials & strcmp(data.participant, subjs{1});
[t1, t2] = meshgrid(find(which_trials_per_subj), find(which_trials_per_subj));

% for each ROI,
% for each run of each subject, get the average RDM
%
for run = 1:metadata.runsPerSubject
    t1_mask = data.runId(t1) == run & data.trialId(t1);
    t2_mask = data.runId(t2) == run & data.trialId(t2);
    run_mask = t1_mask & t2_mask & t1 > t2;
    run_mask_control = repmat(run_mask, 1, 1, numel(control_model_idxs));
    
    for subj = 1:metadata.N
        fprintf('subject %d, run %d\n', subj, run);
        
        for neural_idx = 1:numel(Neural)
            % For a given ROI,
            % For each run for each subject,
            % use a binary mask that says which pairs of trials from the RDM
            % to look at.
            %
            neural_RDM = Neural(neural_idx).RDMs(:,:,subj);
            model_RDM = Model(model_idx).RDMs(:,:,subj);
            control_RDMs = nan(size(neural_RDM, 1), size(model_RDM, 2), numel(control_model_idxs));
            for i = 1:numel(control_model_idxs) 
                control_RDMs(:,:,i) = Model(control_model_idxs(i)).RDMs(:,:,subj);
            end
            
            neural_subRDM = neural_RDM(run_mask);
            model_subRDM = model_RDM(run_mask);
            control_subRDMs = control_RDMs(run_mask_control);
            
            control_subRDMs = reshape(control_subRDMs, size(neural_subRDM, 1), numel(control_model_idxs));
            %rho = partialcorr(neural_subRDM, model_subRDM, control_subRDMs, 'type', 'Spearman');
            rho = corr(neural_subRDM, model_subRDM, 'type', 'Spearman');
            assert(~isnan(rho));

            rhos(run, subj, neural_idx) = rho;
            
            % TODO is this legit? Fisher z-transform?
            %avg_dists(run, subj , neural_idx) = mean(neural_subRDM);
        end
    end
end

%% Get the test choice likelihoods
%
test_log_liks = nan(metadata.runsPerSubject, metadata.N);

% for each run of each subject, get their test choice log likelihood
%
for run = 1:metadata.runsPerSubject    
    for subj = 1:metadata.N
        subj_run_test_trials = data.which_rows & data.runId == run & data.isTest & ~data.timeout &...
            strcmp(data.participant, metadata.allSubjects{goodSubjects(subj)});
        assert(sum(subj_run_test_trials) <= 4 && sum(subj_run_test_trials) >= 3);
        
        
        % OPTION 1 -- do it like we do it in GLM 124 and the KL structure
        % learning effect (kl_structure_learning.m) i.e. take the
        % likelihood for the causal structure corresponding to the current
        % condition only
        %
        model_test_valuess = simulated.valuess(subj_run_test_trials, :); % one for each structure
        model_test_choice_probs = predict(model_test_valuess);
        subj_test_choices = strcmp(data.response.keys(subj_run_test_trials), 'left');
        
        liks = nan(size(model_test_choice_probs));
        for i = 1:size(model_test_choice_probs, 1) % for each test trial
            for j = 1:size(model_test_choice_probs, 2) % for each causal structure
                liks(i,j) = binopdf(subj_test_choices(i), 1, model_test_choice_probs(i, j));
            end
        end
        
        % take average log likelihood !!!
        % this is to account for missing trials where subject timed
        % out. Note that taking the sum would be wrong -- imagine
        % if subject responded on only 1 trial
        test_log_lik = mean(log(liks), 1);
        test_log_lik = test_log_lik(:, 1:3); % exclude M4
        
        M = -1;
        condition = data.condition(subj_run_test_trials);
        condition = condition{1};
        if strcmp(condition, 'irrelevant')
            M = 1;
        elseif strcmp(condition, 'modulatory')
            M = 2;
        else
            assert(strcmp(condition, 'additive'));
            M = 3;
        end
        
        test_log_liks(run, subj) = test_log_lik(M);
        
        % OPTION 2 -- take the overall likelihood; it's basically the same
        % thing since by the end of training, the overall likelihood = the
        % likelihood of the "correct" causal structure
        %
        %{
        model_test_choice_probs = simulated.pred(subj_run_test_trials);
        subj_test_choices = strcmp(data.response.keys(subj_run_test_trials), 'left');
        test_liks = binopdf(subj_test_choices, 1, model_test_choice_probs); % WARNING: won't work on NCF...
        test_log_lik = mean(log(test_liks)); % take the mean -- takes care of TIMEOUTs (think: what if the subject responded on only 1 test trial? we can't just sum them)
        test_log_liks(run, subj) = test_log_lik;
        %test_log_liks(run, subj) = sum(log(test_liks)); % TODO which one is the "right" one?
        %test_log_liks(run, subj) = prod(test_liks);
        %}
    end
end
        
% look at whether this similarity correlates with test behavior
%
r_means = [];
r_sems = [];
all_rs = nan(numel(Neural), metadata.N);

for neural_idx = 1:numel(Neural) % for each ROI
    rhos_roi = rhos(:,:,neural_idx);
    rs = [];
    for subj = 1:metadata.N % for each subject
        x = test_log_liks(:, subj);
        y = rhos_roi(:, subj);
        
        x = zscore(x); % TODO try w/o
        y = zscore(y);
        
        [r, p] = corrcoef(x, y);
        r = r(1,2);
        p = p(1,2);
        fprintf('       ROI = %s, subj = %d, r = %f, p = %f\n', Neural(neural_idx).name, subj, r, p);
        rs = [rs, r];
        all_rs(neural_idx, subj) = r;
        
    end
    
    r_means = [r_means; mean(rs) 0];
    r_sems = [r_sems; sem(rs) 0];        
end

% group-level analysis
%

% For the within-subject analysis, 
% to get the group-level stats you should Fisher z-transform the correlation coefficients 
% (atanh function in matlab) and then do a one-sample t-test against 0.
%
fisher_all_rs = atanh(all_rs);
[h, ps, ci, stats] = ttest(fisher_all_rs');
assert(numel(ps) == numel(Neural));

disp(h);
disp(ps);
n = {Neural.name};
disp(n(logical(h))');
disp(ps(logical(h))');


%{
% plot stuff
%
figure;

subplot(2, 1, 1);
barweb(r_means, r_sems);
ylabel('average r across subjects');
xticklabels({Neural.name});
xtickangle(60);
xlabel('ROI');
%set(gca, 'XTick', []);
title('KL betas correlated with structure learning effect: within-subject analysis', 'Interpreter', 'none');
%}


%% plot
%

figure;

r_means = r_means(:,1);
r_sems = r_sems(:, 1);

h = bar(r_means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
xs = h(1).XData;
hold on;
errorbar(xs, r_means, r_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;

% Put the p-values
%
%bonferroni = false;
for i = 1:numel(Neural)
    significance = @(p) repmat('*', 1, floor(-log10(p)) - 1);
    %{
    if bonferroni
        stars = significance(ps(i) * numel(table_P)); % Bonferroni correction
    else
        stars = significance(ps(i));
    end
    %}
    
    %stars = int2str(length(stars));
    %text(xs(i) - length(stars) * 0.3 / 2, means(i) + sems(i) * 1.2, stars, 'Rotation', 45);
    txt = 'n.s.';
    if ps(i) <= 0.05
        txt = sprintf('p = %.4f', ps(i));
    end
    text(xs(i), r_means(i) + r_sems(i) * 1.2, txt, 'Rotation', 90);
end

% Put the ROI / model names and figure title
%
set(gca, 'xtick', xs);
n = {Neural.name};
labels = {};
for j = 1:numel(n)
    labels{j} = n{j};
end
xticklabels(labels);
xtickangle(55);

ylabel('Avg Fisher z-transformed r');
