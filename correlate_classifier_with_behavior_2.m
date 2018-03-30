% link classifier with behavior: for each ROI, for each subject, plot average per-trial classification accuracy over time -> is it going up over the course of the run?
% Answer: it isn't...
%

which_structures = logical([1 1 0 1 0]);
[data, metadata, simulated, params] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, which_structures);


EXPT = context_expt();
glm = 171;
contrast = 'KL_structures'; 

%EXPT = 'rdms/M1M2M1_4mm/searchlight_tmap_posterior_feedback_onset.nii';
%glm = 0;
%contrast = 'rdms';

%EXPT = 'might/lda_shrinkage_accuracy_countmap_feedback_onset_folds=3_r=2.6667_z-none_use_nosmooth=1_use_tmaps=0.nii';
%glm = 0;
%contrast = 'countmap';


%event = 'trial_onset';
event = 'feedback_onset';
r = 2.6667;

%method = 'cvpatternnet'; % NN; BEST performance on the peak LDA voxels; not so much on the GNB peak voxels => LDA is indeed better
method = 'cvfitcnb'; % gaussian naive bayes, just like searchmight's gnb_searchmight; MEDIUM performance
%method = 'cvfitcdiscr'; % linear discriminant analysis, similar to searchmight's lda_shrinkage but not quite; WORST performance
runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none';

use_tmaps = false;
use_nosmooth = true; 

%classifier = 'gnb_searchmight';
classifier = 'lda_shrinkage';

p = 0.001;
direct = '+';
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3


n_iter = 3; % how many iterations for each subject
trial_cutoff = 1; % average classifier posterior from trials >= trial_cutoff in each run (to smoothen the noise) 
best_k_voxels = 3; % look at best k voxels in each ROI

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

accs = [];

clear cached_o;
cached_filename = fullfile('temp', sprintf('corr_class_w_behav_%s_niter=%d_bestk=%d_contrast=%s_event=%s.mat', method, n_iter, best_k_voxels, contrast, event));

% load cached outputs, to avoid having to re-classify
% WARNING: make sure it's the correct contrast, etc; basically, never assume you've cached the right thing
%
load_cached = false;
if load_cached
    load(cached_filename, 'cached_o');
end


for i = 1:size(region, 1) % for each ROI
%for i = 1:6 % for each ROI
    fprintf('ROI = %s\n', region{i});

    clust_idx = CI(cor(i,1), cor(i,2), cor(i,3));
    clust_mask = CI == clust_idx;
    clust_vox = find(clust_mask);

    %clust_mask = mask; HACK to do the whole-brain
    %clust_vox = find(mask);

    liks_classifier = []; % avg test log lik for each run for each subject, according to classifier posterior
    liks_model = []; % same but according to model posterior

    acc = [];
    subj_idx = 0;
    for subj = subjs % for each subject
        subject = metadata.allSubjects(subj); % 'con001' ... 'con025'
        subj_idx = subj_idx + 1;

        if ~load_cached % optionally, the outputs are precomputed and cached
            filename = sprintf(file_format, classifier, event, subj, r, z_score, use_nosmooth, use_tmaps);
            [~, ~, amap] = load_mask(filename); % get accuracy map for subject

            [~, j] = sort(amap(clust_mask), 'descend');
            j = j(1:best_k_voxels); % get k voxels with peak accuracy within the ROI
            [x,y,z] = ind2sub(size(clust_mask), clust_vox(j));
            fprintf('\nsubj = %d, max acc(s) = %s\n', subj, sprintf('%.2f, ', amap(sub2ind(size(amap),x,y,z))));

            % average "votes" across top k voxels in ROI
            %
            o = zeros(size(targets));
            for k = 1:best_k_voxels % for each of the top k voxels
                % get sphere around voxel for the given ROI (same sphere as the searchmight)
                %
                sphere_mask = create_spherical_mask_helper(mask, x(k), y(k), z(k), r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);
                sphere_activations = get_activations_submask(sphere_mask, activations);

                % convert to inputs/targets
                %
                [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, [subj], sphere_activations, predict_what, z_score, data, metadata);
                inputs = inputs(:, sum(inputs,1) ~= 0); % clear columns of all 0's, previously NaNs b/c of use_nosmooth; we clear those out in might_getNeighbors() for searchmight

                % run classifier n_iter times, to account for randomness in the folds. Average outputs (posteriors) 
                %
                for iter = 1:n_iter
                    [~, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, [subj], []);
                    o = o + outputs;
                end
            end

            o = o / (n_iter * best_k_voxels); % final output for ROI (for that subject) = average across iterations of top k voxels
            cached_o{i}{subj} = o;
        else
            % load cached outputs
            o = cached_o{i}{subj};
        end

        accuracy = classify_get_accuracy(o, targets); % "average" accuracy, according to average votes across iterations

        %if accuracy < 35
        %    fprintf('SKIPPING SUBJECT %d: accuracy too low (%.2f)\n', subj, accuracy);
        %    continue;
        %end

        rid = data.runId(which_rows);
        tid = data.trialId(which_rows);
        subplot(size(region, 1), numel(subjs), (i-1) * numel(subjs) + subj_idx); 
        a = [];
        e = [];
        otar = o .* targets; 
        for t = trials
            a_ts = sum(otar(tid == t, :), 2);
            a_t = mean(a_ts);
            a = [a a_t];
            e = [e std(a_ts) / sqrt(numel(a_ts))];
        end
        errorbar(a, e);


        % use outputs as posteriors in simulations, compare w/ model posterior on test log lik
        %
        %{
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
            P_n = zeros(size(which_structures));

            % average classifier posterior from last X trials ~= structure posterior
            P_n(which_structures) = mean(o(rid == run & tid >= trial_cutoff, :), 1); % the important part -- notice we use the average outputs of the n_iter iterations

            % option #1: this is a very noisy estimate; use it only to bias the model posterior (as opposed to replace it)
            P_n_model = simulated.P(which_run & data.trialId == metadata.trainingTrialsPerRun, :);
            P_n = P_n_model + P_n * 0.5; % TODO param, free?
            P_n = P_n / sum(P_n);

            % option #2: set P_n = MAP structure
            %[~,map] = max(P_n);
            %P_n = zeros(size(P_n));
            %P_n(map) = 1;

            train_results.P_n = P_n;
            train_results.ww_n = ww_n;
            test_results = model_test(test_x, test_k, train_results, params); % simulate run using P_n ~= posterior from classifier

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
        %}

        acc = [acc, accuracy];
    end

    %{
    [h, p, ci, stat] = ttest(liks_classifier, liks_model) % paired sample t-test TODO transform log likelihoods to make them more Gaussian-y
    ttest_ts = [ttest_ts; stat.tstat];
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(liks_classifier) mean(liks_model)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
    %}

    accs = [accs; acc];
end


if ~load_cached
    save(cached_filename); % save cached outputs, to avoid having to reclassify each time
end
