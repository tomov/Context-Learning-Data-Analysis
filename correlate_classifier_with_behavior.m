% using the classifier posterior from various ROIs as proxy for structure posterior; see if it's better at predicting test choices than model posterior
%
% how to do it -- 
% 0) for each ROI (whole cluster) in the KL_structures GLM (maybe RSA too?)
% 1) for each subject, pick highest accuracy voxel in that cluster
% 2) use the sphere around that voxel (i.e. what got that accuracy in the first place) 
% 3) run your own classifier (to get the posterior over conditions at each time point; searchmight didn't give you that)
% 4) predict choices based on that as posterior (log lik) 
% 6) see if that's better, on average, than the model log lik


which_structures = logical([1 1 0 1 0]);
[data, metadata, simulated, params] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, which_structures);


%EXPT = context_expt();
EXPT = 'rdms/M1M2M1_4mm/searchlight_tmap_posterior_feedback_onset.nii';
%glm = 171;
glm = 0;
%contrast = 'KL_structures'; 
contrast = 'rdms';

event = 'trial_onset';
r = 2.6667;

method = 'cvfitcnb'; % gaussian naive bayes, just like searchmight 
runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none';

n_iter = 10; % how many iterations for each subject

use_tmaps = false; % <-- slightly better if true; but stick with betas for consistency w/ RDMs
use_nosmooth = true; 

%classifier = 'gnb_searchmight';
classifier = 'lda_shrinkage';

p = 0.001;
direct = '+';
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end

[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glm, contrast, p, direct, alpha, Dis, Num);

file_format = 'might/%s_accuracy_%s_subj=%d_folds=3_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii';
maskfile = 'masks/mask.nii';

[mask, Vmask] = load_mask(maskfile);
Vmask.fname = 'temp.nii'; % change immediately! in case we screw it up
[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);

%activations = get_activations(maskfile, event, data, metadata, use_nosmooth);

ttest_means = [];
ttest_sems = [];
ttest_ps = [];
ttest_ts = [];

for i = 1:size(region, 1) % for each ROI
    fprintf('ROI = %s\n', region{i});

    clust_idx = CI(cor(i,1), cor(i,2), cor(i,3));
    clust_mask = CI == clust_idx;
    clust_vox = find(clust_mask);

    %clust_mask = mask; HACK to do the whole-brain
    %clust_vox = find(mask);

    liks_classifier = []; % avg test log lik for each run for each subject, according to classifier posterior
    liks_model = []; % same but according to model posterior

    acc = [];
    for subj = subjs % for each subject
        subject = metadata.allSubjects(subj); % 'con001' ... 'con025'

        filename = sprintf(file_format, classifier, event, subj, r, z_score, use_nosmooth, use_tmaps);
        [~, ~, amap] = load_mask(filename); % get accuracy map

        [~, j] = max(amap(clust_mask));
        [x,y,z] = ind2sub(size(clust_mask), clust_vox(j));
        fprintf('\nsubj = %d, max acc = %.4f\n', subj, amap(x,y,z));

        % get sphere around peak accuracy voxel for the given cluster (same sphere as the searchmight)
        %
        sphere_mask = create_spherical_mask_helper(mask, x, y, z, r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);
        sphere_activations = get_activations_submask(sphere_mask, activations);

        % convert to inputs/targets
        %
        [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, [subj], sphere_activations, predict_what, z_score, data, metadata);
        inputs = inputs(:, sum(inputs,1) ~= 0); % clear columns of all 0's, previously NaNs b/c of use_nosmooth; we clear those out in might_getNeighbors() for searchmight


        % run classifier n_iter times, to account for randomness in the folds. Average outputs (posteriors) 
        %
        o = zeros(size(targets));
        for iter = 1:n_iter
            [~, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, [subj], []);
            o = o + outputs;
        end
        o = o / n_iter;
        accuracy = classify_get_accuracy(o, targets);

        %if accuracy < 35
        %    fprintf('SKIPPING SUBJECT %d: accuracy too low (%.2f)\n', subj, accuracy);
        %    continue;
        %end

        % use outputs as posteriors in simulations, compare w/ model posterior on test log lik
        %
        rid = data.runId(which_rows);
        tid = data.trialId(which_rows);
        for run = runs % for each run
            which_run = data.which_rows & strcmp(data.participant, subject) & data.runId == run;
            assert(sum(which_run) == metadata.trialsPerRun);

            which_run_test = data.which_rows & ~data.isTrain & ~data.timeout & strcmp(data.participant, subject) & data.runId == run;

            % use the classifier output as posterior to compute the test log likelihood
            %
            [train_x, train_k, train_r, test_x, test_k] = convert_run(data, metadata, subject, run);
            ww_n = simulated.ww_n{which_run & data.trialId == metadata.trainingTrialsPerRun}; % notice we still use the weights from the model
            %P_n = simulated.P(which_run & data.trialId == metadata.trainingTrialsPerRun,:); <-- sanity; pred_classifier should be same as pred_model in that case
            P_n = mean(outputs(rid == run & tid > 0, :), 1); % the important part
            train_results.P_n = P_n;
            train_results.ww_n = ww_n;

            test_results = model_test(test_x, test_k, train_results, params);

            pred_classifier = test_results.choices(which_run_test(data.which_rows & ~data.isTrain & strcmp(data.participant, subject) & data.runId == run)); % P(choose sick) on the test trials, according to the classifier


            % compare with the test log likelihood according to the model
            %
            pred_model = simulated.pred(which_run_test); % P(choose sick) on the test trials, according to model
            X = data.chose_sick(which_run_test); % actual subject choices

            test_liks_classifier = binopdf(X, 1, pred_classifier); 
            test_liks_model = binopdf(X, 1, pred_model); 
            assert(numel(test_liks_classifier) == numel(X));
            assert(numel(test_liks_model) == numel(X));

            avg_loglik_classifier = mean(log(test_liks_classifier));
            avg_loglik_model = mean(log(test_liks_model));

            liks_classifier = [liks_classifier, avg_loglik_classifier];
            liks_model = [liks_model, avg_loglik_model];
        end

        acc = [acc, accuracy];
    end

    [h, p, ci, stat] = ttest(liks_classifier, liks_model) % paired sample t-test TODO transform log likelihoods to make them more Gaussian-y
    ttest_ts = [ttest_ts; stat.tstat];
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(liks_classifier) mean(liks_model)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
end
