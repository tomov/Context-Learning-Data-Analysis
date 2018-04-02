glmodel = 163;
regressor = 'KL_structures';
contrast = 'KL_structures - KL_weights';
which_structures = logical([1 1 0 1 0]);
[data, metadata, simulated, params] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, which_structures);

%function [means, sems, ps, ts] = betas_to_behavior_2(glmodel, regressor, params, which_structures, contrast)
% Correlate peak voxels from different contrasts with behavior across subjects (akin to correlate_classifier_with_behavior_3.m).
% Picks max beta in each subject (avg across runs)
%
% INPUT:
% glmodel = glm as in context_create_multi.m
% regressor = name of regressor to get betas for
% [params, which_structures] = model_default_params();
% contrast = optional contrast from which to extract the clusters; by
%            default set to regressor 
%
% EXAMPLES:
% [params, which_structures] = model_params('results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat')
% [means, sems, ps, ts] = betas_to_behavior(163, 'KL_structures', 'sphere', params, which_structures, 'KL_structures - KL_weights')
%

if ~exist('contrast', 'var') || isempty(contrast)
    contrast = regressor;
end

EXPT = context_expt();
%glmodel = 154;
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
%r = 1.814;
r = 2.6667;
direct = '+';
%regressor = 'KL_structures'; % contrasts not supported yet b/c ot load_run_betas
%what = 'voxel';

best_k_voxels = 1;



%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Find the peak voxels 
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);


rhos = [];
ps = [];

figure;

% Average betas from the cluster of each peak voxel
% see correlate_classifier_with_behavior_3.m
%
for i = 1:numel(region)
    fprintf('\nGetting betas for cluster t=%.3f extent=%d roi=%s peak=[%d %d %d]\n', stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));

    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);

    % create mask from cluster
    clust_idx = CI(x,y,z);
    clust_mask = CI == clust_idx;

    % take voxels from cluster and convert to MNI
    clust_vox = find(clust_mask);
    assert(numel(clust_vox) == extent(i));
    [x, y, z] = ind2sub(size(clust_mask), clust_vox);
    voxels = [x, y, z];

    %betas = load_run_betas(glmodel, regressor, mni(i,:));
    %betas = mean(betas, 2);

    % get the corresponding betas
    clust_betas = load_run_betas(glmodel, regressor, voxels, false); % n_subjects x n_runs x n_voxels
    % average betas across runs
    clust_betas = squeeze(mean(clust_betas, 2)); % n_subjects x n_voxels

    subj_betas = [];
    subj_logliks = [];
    for subj_idx = 1:metadata.N % for each subject
        subject = metadata.subjects(subj_idx); % 'con001' ... 'con025'

        [~, j] = sort(clust_betas(subj_idx, :), 'descend');
        j = j(1:best_k_voxels); % get k voxels with max betas within the ROI
        vox_betas = clust_betas(subj_idx, j); % betas of top k voxels in that ROI for that subject

        %vox_betas = betas(subj_idx);

        run_logliks = [];
        for run = 1:metadata.runsPerSubject % for each run, compute test log likelihood
            which_run_test = data.which_rows & ~data.isTrain & ~data.timeout & strcmp(data.participant, subject) & data.runId == run;

            pred = simulated.pred(which_run_test); % P(choose sick) on the test trials, according to model
            X = data.chose_sick(which_run_test); % actual subject choices
            test_liks = binopdf(X, 1, pred); 
            assert(numel(test_liks) == numel(X));
            avg_loglik = mean(log(test_liks));

            run_logliks = [run_logliks, avg_loglik];
        end

        subj_betas = [subj_betas, mean(vox_betas)];
        subj_logliks = [subj_logliks, mean(run_logliks)];
    end


    [rho, p] = corr(subj_betas', subj_logliks', 'type', 'Pearson');
    rhos = [rhos, rho];
    ps = [ps, p];

    fprintf('r = %.4f, p = %.4f\n', rho, p);

    subplot(1, numel(region), i);
    scatter(subj_betas, subj_logliks);
    lsline;
    xlabel('beta');
    ylabel('loglik');
    title(region{i}, 'interpreter', 'none');

end
        


