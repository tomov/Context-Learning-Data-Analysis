function [means, sems, ps] = betas_to_behavior(glmodel, regressor, what)

% Correlate peak voxels from different contrasts with behavior
% Also tries with clusters of voxels (average the betas)
%
% EXAMPLES:
% betas_to_behavior(123, 'surprise', 'voxel')
% betas_to_behavior(148, 'KL_weights', 'voxel')
% betas_to_behavior(154, 'KL_structures', 'voxel')
% betas_to_behavior(154, 'KL_weights', 'voxel')
% betas_to_behavior(154, 'KL_structures', 'sphere')
% betas_to_behavior(154, 'KL_weights', 'sphere')
% betas_to_behavior(154, 'KL_structures', 'cluster')
% betas_to_behavior(154, 'KL_weights', 'cluster')
%

EXPT = context_expt();
%glmodel = 154;
p = 0.001;
alpha = 0.001;
Dis = 20;
Num = 1;
r = 1.814;
direct = '+';
%regressor = 'KL_structures'; % contrasts not supported yet b/c ot load_run_betas
%what = 'voxel';

assert(ismember(what, {'voxel', 'sphere', 'cluster'}));

filename = ['results/betas_to_behavior_glm', num2str(glmodel), '_', regressor, '_', what, '.mat'];


%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Find the peak voxels 
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glmodel, regressor, p, direct, alpha, Dis, Num);


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
[test_liks, test_RTs] = get_test_behavior();

save(filename);
fprintf('SAVING to %s\n', filename);


%% within-subject analysis using a linear mixed effects model and/or t-tests
%
load(filename);

[means, sems, ps] = correlate_neural_and_behavior(betas, region, test_liks, [regressor, ' ', what, ' betas from GLM ', num2str(glmodel), ' correlated with test log likelihood: t-test']);
