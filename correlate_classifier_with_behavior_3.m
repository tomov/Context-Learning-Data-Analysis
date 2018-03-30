% For each ROI, see if average accuracy correlates with test log likelihood (across subjects)
%

which_structures = logical([1 1 0 1 0]);
[data, metadata, simulated, params] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, which_structures);


%EXPT = context_expt();
%glm = 171;
%contrast = 'KL_structures'; 

EXPT = 'rdms/M1M2M1_4mm/searchlight_tmap_posterior_feedback_onset.nii';
glm = 0;
contrast = 'rdms';

%EXPT = 'rdms/M1M2M1_4mm/searchlight_tmap_prior_trial_onset.nii';
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
best_k_voxels = 1; % look at best k voxels in each ROI

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


rhos = [];
ps = [];

for i = 1:size(region, 1) % for each ROI
    fprintf('\nROI = %s\n', region{i});

    clust_idx = CI(cor(i,1), cor(i,2), cor(i,3));

    %clust_mask = mask; HACK to do the whole-brain
    %clust_vox = find(mask);

    liks_classifier = []; % avg test log lik for each run for each subject, according to classifier posterior
    liks_model = []; % same but according to model posterior

    subj_accs = [];
    subj_logliks = [];
    subj_idx = 0;
    for subj = subjs % for each subject
        subject = metadata.allSubjects(subj); % 'con001' ... 'con025'
        subj_idx = subj_idx + 1;

        filename = sprintf(file_format, classifier, event, subj, r, z_score, use_nosmooth, use_tmaps);
        [~, ~, amap] = load_mask(filename); % get accuracy map for subject

        clust_mask = CI == clust_idx & ~isnan(amap);
        clust_vox = find(clust_mask);

        [~, j] = sort(amap(clust_mask), 'descend');
        j = j(1:best_k_voxels); % get k voxels with peak accuracy within the ROI
        [x,y,z] = ind2sub(size(clust_mask), clust_vox(j));
        vox_accs = amap(sub2ind(size(amap),x,y,z));
        %fprintf('subj = %d, max acc(s) = %s\n', subj, sprintf('%.2f, ', vox_accs));

        rid = data.runId(which_rows);
        tid = data.trialId(which_rows);
        run_logliks = [];
        for run = runs % for each run, compute test log likelihood
            which_run_test = data.which_rows & ~data.isTrain & ~data.timeout & strcmp(data.participant, subject) & data.runId == run;

            pred = simulated.pred(which_run_test); % P(choose sick) on the test trials, according to model
            X = data.chose_sick(which_run_test); % actual subject choices
            test_liks = binopdf(X, 1, pred); 
            assert(numel(test_liks) == numel(X));
            avg_loglik = mean(log(test_liks));

            run_logliks = [run_logliks, avg_loglik];
        end

        subj_accs = [subj_accs, mean(vox_accs)];
        subj_logliks = [subj_logliks, mean(run_logliks)];
    end

    [rho, p] = corr(subj_accs', subj_logliks', 'type', 'Pearson');
    rhos = [rhos, rho];
    ps = [ps, p];

    fprintf('r = %.4f, p = %.4f\n', rho, p);
end

