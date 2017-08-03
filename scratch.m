


[cluster,Vcluster] = load_mask('masks/glm0_light_cluster_t=5.435_extent=17_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii');
[edited,Vedited] = load_mask('masks/light_prior_trial-onset_L_dlPFC1_ClusterMask_searchlight_tmap_x=-36_y=10_z=24_17voxels_edited.nii');

[x,y,z] = ind2sub(size(edited), find(edited));
Cedited = [x y z];

[x,y,z] = ind2sub(size(cluster), find(cluster));
Ccluster = [x y z];

MNIedited = cor2mni(Cedited, Vedited.mat);
MNIcluster = cor2mni(Ccluster, Vcluster.mat);

MNIedited = sort(MNIedited)
MNIcluster = sort(MNIcluster)

assert(isequal(MNIcluster, MNIedited));


%betas1 = get_betas('masks/glm0_light_sphere_t=5.435_extent=24_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii', 'feedback_onset', data, metadata, false);


%{
clear all;
close all;
[means{1}, sems{1}, ps{1}] = betas_to_behavior(123, 'surprise', 'voxel');
[means{2}, sems{2}, ps{2}] = betas_to_behavior(148, 'KL_weights', 'voxel');
[means{3}, sems{3}, ps{3}] = betas_to_behavior(154, 'KL_structures', 'voxel');
[means{4}, sems{4}, ps{4}] = betas_to_behavior(154, 'KL_weights', 'voxel');
[means{5}, sems{5}, ps{5}] = betas_to_behavior(154, 'KL_structures', 'sphere');
[means{6}, sems{6}, ps{6}] = betas_to_behavior(154, 'KL_weights', 'sphere');
[means{7}, sems{7}, ps{7}] = betas_to_behavior(154, 'KL_structures', 'cluster');
[means{8}, sems{8}, ps{8}] = betas_to_behavior(154, 'KL_weights', 'cluster');

save('results/betas_to_behavior_alpha=0.001_Num=1.mat');

%}

%{


% eLife reviewer comment -- contextual tracking confounding KL

%% Load data
%
[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

context_changed = data.contextId ~= circshift(data.contextId, 1);
context_changed(data.trialId == 1 & data.isTrain) = 0; % TODO 0 or 1?

cue_changed = data.cueId ~= circshift(data.cueId, 1);
cue_changed(data.trialId == 1 & data.isTrain) = 0; % TODO 0 or 1?

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


simulated = simulate_subjects(data, metadata, params, which_structures);

%% mad stats
%
which = which_trials & data.trialId > 1; % ignore trial 1

KL_c_change = simulated.surprise(which & context_changed);
KL_no_c_change = simulated.surprise(which & ~context_changed);

KL_x_change = simulated.surprise(which & cue_changed);
KL_no_x_change = simulated.surprise(which & ~cue_changed);

% are they from same distribution?
[h, p, ks2stat] = kstest2(KL_c_change, KL_no_c_change);
fprintf('Kolmogorov-Smirnov for context changes: p(same distribution) = %e\n', p);

[h, p, ks2stat] = kstest2(KL_x_change, KL_no_x_change);
fprintf('Kolmogorov-Smirnov for cue changes: p(same distribution) = %e\n', p);

% do higher KL's correlate with changes in context?
[rho, p] = corr([simulated.surprise(which), context_changed(which)], 'type', 'Spearman');
rho = rho(1,2);
p = p(1,2);
fprintf('Spearman for context changes: rho = %f, p(no correlation) = %e\n', rho, p);

[rho, p] = corr([simulated.surprise(which), cue_changed(which)], 'type', 'Spearman');
rho = rho(1,2);
p = p(1,2);
fprintf('Spearman for cue changes: rho = %f, p(no correlation) = %e\n', rho, p);

% do KL's for context change vs. no context change have the same means?
% (wrong: assumes they're gaussian; they're not...)
[p, anovatab, stats] = anova1(simulated.surprise(which), context_changed(which));
fprintf('ANOVA for context changes: p(means are equal) = %e\n', p);

[p, anovatab, stats] = anova1(simulated.surprise(which), cue_changed(which));
fprintf('ANOVA for cue changes: p(means are equal) = %e\n', p);

%}



%{
batch_size = 1000;
r = 1.814;
for batch=1:217
    fprintf('Batch = %d\n', batch);
    start_idx = (batch - 1) * batch_size + 1;
    end_idx = batch * batch_size;
    rdms_searchlight(start_idx, end_idx, 1.814);
end
%}


%{

method = 'cvglmnet';
mask = fullfile('masks', 'hippocampus.nii');
z_score = 'z-none';
predict_what = 'condition';
runs = [1:9];
trials = [1:24];
subjs = getGoodSubjects();

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, subjs);

leftout_run = 9;
runs(runs == leftout_run) = [];
which_rows = data.which_rows & ismember(data.trialId, trials) & data.runId ~= leftout_run;

m = containers.Map(metadata.subjects, subjs);
subj_id = cellfun(@(x) m(x), data.participant(which_rows));

foldid = data.runId(which_rows) + subj_id * (metadata.runsPerSubject + 1);
assert(numel(unique(foldid)) == numel(subjs) * (metadata.runsPerSubject - 1));
unique_foldid = unique(foldid);
foldid = arrayfun(@(x) find(unique_foldid == x), foldid);
assert(numel(unique(foldid)) == numel(subjs) * (metadata.runsPerSubject - 1));
assert(max(foldid) == numel(subjs) * (metadata.runsPerSubject - 1));


[classifier, inputs, targets, outputs, which_rows] = classify_train(method, runs, trials, subjs, mask, predict_what, z_score, foldid);

[test_inputs, test_targets, test_outputs, test_which_rows] = classify_test(method, classifier, leftout_run, trials, subjs, mask, predict_what, z_score);

%}

%{
EXPT = context_expt();
subj = 1;
model = 139;
modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
load(fullfile(modeldir,'SPM.mat'));

[N,K] = size(SPM.xX.X)
%}


%{
o = all_outputs(all_which_rows, :);
t = all_targets(all_which_rows, :);
acc = classify_get_accuracy(ones(size(o)) / 3, t)
%}


%{
handle = figure;
        %set(handle, 'Position', [500, 500, 450, 200])
        
        subplot(2, 1, 1);

        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.isTrain & data.trialId == 20 & strcmp(data.contextRole, condition);
            
            P = simulated.P(which_rows, which_structures);             
            P_means = [P_means; mean(P, 1)];
        end

        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Posterior probability');
        legend({'M1', 'M2', 'M3'}, 'Position', [0.2 0.3 1 1]);
        title('Model');
    
        %
        % Figure 3B: Choice probabilities on test trials for model vs. humans
        %
        
        subplot(2, 1, 2);
        
        % TODO SAM superimpose human dots w/ error bars; rm error bars from model
        %
        % TODO SAM also SEMs -- within-subject errors: wse.m (in the dropbox)
        % TODO SAM also mytitle.m -- easier to plot left-oriented panel titles
        
        % Choice probabilities for model
        %
        human_means = [];
        human_sems = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
            
            x1c1 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 0);
            x1c3 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 2);
            x3c1 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 0);
            x3c3 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 2);

            human_means = [human_means; mean(x1c1) mean(x1c3) mean(x3c1) mean(x3c3)];
            human_sems = [human_sems; sem(x1c1) sem(x1c3) sem(x3c1) sem(x3c3)];
        end

        hold on;
        %errorbar(human_means(1,:), human_sems(1,:));
        hold off;
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Choice probability');
        legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 0.2 1 1]);
        title('Human subjects');
%}
