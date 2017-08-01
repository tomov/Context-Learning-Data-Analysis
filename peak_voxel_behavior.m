% Correlate peak voxels from different contrasts with behavior
%

close all;
clear all;

EXPT = context_expt();
glmodel = 154;

%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Simulate behavior
%
prior_variance = 0.1249;
inv_softmax_temp = 2.0064;
params = [prior_variance inv_softmax_temp];
which_structures = logical([1 1 1 0]);
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows, false);

%% Some random voxels as controls
%
rand_voxels = gen_rand_voxels(fullfile('masks', 'mask.nii'), 20);

%% Find the peak voxels 
%
[KL_structures_voxels, KL_structures_rois] = get_peak_voxels_from_contrast(EXPT, glmodel, 'KL_structures');
[KL_weights_voxels, KL_weights_rois] = get_peak_voxels_from_contrast(EXPT, glmodel, 'KL_weights');

%% Get the betas
%
KL_structures_betas = load_run_betas(glmodel, 'KL_structures', [KL_structures_voxels; rand_voxels]);
KL_weights_betas = load_run_betas(glmodel, 'KL_weights', [KL_weights_voxels; rand_voxels]);

%% Get the behavioral correlates
%
[test_liks, test_RTs] = get_test_behavior();

save('results/peak_voxel_behavior.mat');

%% within-subject analysis using a linear mixed effects model and/or t-tests
%

load('results/peak_voxel_behavior.mat');

correlate_neural_and_behavior(KL_weights_betas, [KL_weights_rois', repmat({'random'}, 1, size(rand_voxels, 1))], test_liks, 'KL_weight betas correlated with test log likelihood: t-test');
correlate_neural_and_behavior(KL_structures_betas, [KL_structures_rois', repmat({'random'}, 1, size(rand_voxels, 1))], test_liks, 'KL_structure betas correlated with test log likelihood: t-test');
