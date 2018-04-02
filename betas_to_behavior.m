function [means, sems, ps, ts] = betas_to_behavior(glmodel, regressor, what, params, which_structures, contrast)

% Correlate peak voxels from different contrasts with behavior within-subjects (i.e. across runs). See if correlation > 0 across subjects.
% Also tries with clusters of voxels (average the betas)
%
% INPUT:
% glmodel = glm as in context_create_multi.m
% regressor = name of regressor to get betas for
% what = 'voxel', 'sphere', or 'cluster' -- what area to take around the
%        peak voxel from each cluster
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

assert(ismember(what, {'voxel', 'sphere', 'cluster'}));

filename = ['results/betas_to_behavior_glm', num2str(glmodel), '_', regressor, '_', what, '_', contrast, '.mat'];


%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Find the peak voxels 
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);


%% Get the betas
%

% For each peak voxel, take the whole cluster around it and average the
% betas
%
betas = nan(metadata.N, metadata.runsPerSubject, numel(region));

switch what
    case 'voxel'
        % Take betas from peak voxels only
        %
        betas = load_run_betas(glmodel, regressor, mni);
        
        
    case 'sphere'
        % Average betas from a sphere around each peak voxel
        % TODO maybe intersect with the cluster too
        %
        for i = 1:numel(region)
            fprintf('\nGetting betas for sphere of peak voxel t=%.3f extent=%d roi=%s peak=[%d %d %d]\n', stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));

            % create spherical mask around voxel
            [~, voxels] = create_spherical_mask(mni(i,1), mni(i,2), mni(i,3), r);

            % get the corresponding betas
            sphere_betas = load_run_betas(glmodel, regressor, voxels);

            % average them
            fprintf('Averaging %d betas\n', size(sphere_betas, 3));
            sphere_betas = mean(sphere_betas, 3);

            betas(:,:,i) = sphere_betas;
        end
        
        
    case 'cluster'
        % Average betas from the cluster of each peak voxel
        %
        for i = 1:numel(region)
            fprintf('\nGetting betas for cluster t=%.3f extent=%d roi=%s peak=[%d %d %d]\n', stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));

            x = cor(i,1);
            y = cor(i,2);
            z = cor(i,3);
            assert(immse(stat(i), Y(x,y,z)) < 1e-6);

            % create mask from cluster
            clust_idx = CI(x,y,z);
            mask = CI == clust_idx;

            % take voxels from cluster and convert to MNI
            voxels = find(mask);
            assert(numel(voxels) == extent(i));
            [x, y, z] = ind2sub(size(mask), voxels);
            voxels = cor2mni([x, y, z], V.mat);

            % get the corresponding betas
            clust_betas = load_run_betas(glmodel, regressor, voxels);

            % average them
            fprintf('Averaging %d betas\n', size(clust_betas, 3));
            clust_betas = mean(clust_betas, 3);

            betas(:,:,i) = clust_betas;
        end
        
    otherwise 
        assert(false);
end


%% Get the behavioral correlates
%
[test_liks, test_RTs] = get_test_behavior(params, which_structures);


%% within-subject analysis using a linear mixed effects model and/or t-tests
%
[means, sems, ps, ts] = correlate_neural_and_behavior(betas, region, test_liks, [regressor, ' ', what, ' betas from GLM ', num2str(glmodel), ' correlated with test log likelihood: t-test']);


save(filename);
fprintf('SAVING to %s\n', filename);
