% Correlate peak voxels from different contrasts with behavior
% Also tries with clusters of voxels (average the betas)
% TODO dedupe with peak_voxel_behavior.m
%

close all;
clear all;

EXPT = context_expt();
glmodel = 154;
p = 0.001;
direct = '+';
regressor = 'KL_structures'; % contrasts not supported yet b/c ot load_run_betas

filename = ['results/betas_to_behavior_glm', num2str(glmodel), '_', regressor];

%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Find the peak voxels 
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glmodel, regressor, p, direct);


%[KL_structures_voxels, KL_structures_rois] = get_peak_voxels_from_contrast(EXPT, glmodel, 'KL_structures');
%[KL_weights_voxels, KL_weights_rois] = get_peak_voxels_from_contrast(EXPT, glmodel, 'KL_weights');

%% Get the betas
%

% For each peak voxel, take the whole cluster around it and average the
% betas
%
betas = nan(metadata.N, metadata.runsPerSubject, numel(region));

for i = 1:numel(region)
    fprintf('Getting betas for cluster t=%.3f extent=%d roi=%s peak=[%d %d %d]\n', stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));
    
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
    clust_betas = mean(clust_betas, 3);
    
    betas(:,:,i) = clust_betas;
end

%KL_structures_betas = load_run_betas(glmodel, 'KL_structures', [KL_structures_voxels; rand_voxels]);
%KL_weights_betas = load_run_betas(glmodel, 'KL_weights', [KL_weights_voxels; rand_voxels]);

%% Get the behavioral correlates
%
[test_liks, test_RTs] = get_test_behavior();

save(filename);
fprintf('SAVING to %s\n', filename);

%% within-subject analysis using a linear mixed effects model and/or t-tests
%

load(filename);

correlate_neural_and_behavior(betas, regions, test_liks, [regressor, ' betas correlated with test log likelihood: t-test']);
