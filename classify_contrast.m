function [rois, targets, which_rows] = classify_contrast(EXPT, glmodel, contrast, what, event)

% Train classifier using ROIs from a contrast.
% similar to classify_searchlight.m
%
% INPUT:
% EXPT = expt struct
% glmodel = glm as in context_create_multi.m
% contrast = contrast from which to extract the clusters; by
% what = 'sphere', or 'cluster' -- what area to take around the
%        peak voxel from each cluster
% EXAMPLE:
% classify_contrast(context_expt(), 171, 'KL_structures', 'sphere')
%
% OUTPUT:
% rois = struct array of results for each ROI
%

use_tmaps = false;
use_nosmooth = true;

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end

dirname = 'classify';

% fmri cluster params
%
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
%r = 1.814;
r = 2.6667; 
%r = 6.6667; 
direct = '+';

% classifier params
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
%method = 'cvglmnet'; <-- sucks on synthetic data
method = 'cvpatternnet';
runs = 1:metadata.runsPerSubject;
trials = 6:metadata.trainingTrialsPerRun; % TODO arbitrary
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none';
%z_score = 'z-run';
%z_score = 'z-run-voxel';
%z_score = 'pca-subj';
%z_score = 'z-run-voxel-pca-subj';
%event = 'feedback_onset';

% get whole-brain betas
%
whole_brain_activations = get_activations('masks/mask.nii', event, data, metadata, use_nosmooth);

% extract clusters
%
switch what
    case 'sphere' 
        [filenames, masknames] = create_sphere_masks_from_contrast(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, r);

    case 'cluster'
        [filenames, masknames] = create_masks_from_contrast(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);

    otherwise 
        assert(false);
end


% train & cross-validate classifiers
%
for mask_idx = 1:numel(filenames)
    fprintf('Cluster %d: mask %s\n', mask_idx, filenames{mask_idx});

    % get betas
    %
    [m, V] = load_mask(filenames{mask_idx});
    activations = get_activations_submask(m, whole_brain_activations);

    % if a voxel is NaN even for 1 trial, ignore it
    %
    good_voxels = sum(isnan(activations(data.which_rows & data.isTrain,:)), 1) == 0; 
    if sum(good_voxels) == 0
        assert(use_nosmooth); % doesn't happen if we have smoothing b/c we are using the mask that by definition contains no NaNs
        warning('Skipping cluster -- no good voxels');
        continue;
    end
    if sum(good_voxels) < size(activations, 2)
        assert(use_nosmooth);
        activations = activations(:,good_voxels);
        warning(sprintf('Cluster has only %d good voxels', sum(good_voxels)));
    end

    % one classifier for all subjects -- not standard
    %
    %[inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs, activations, predict_what, z_score, data, metadata);
    %[classifier, outputs, accuracy] = classify_train_helper(method, inputs, targets, runs, trials, subjs, []);

    %[classifier, ~, targets, outputs, which_rows, accuracy] = classify_train(method, runs, trials, subjs, filenames{mask_idx}, predict_what, z_score, event); <-- too much overhead loading the betas

    % one classifier per subject
    %
    accuracies = [];
    ps = [];
    for subj = subjs
        fprintf('    subj %d\n', subj);

        [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, [subj], activations, predict_what, z_score, data, metadata);
        [classifier, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, [subj], []);
        accuracies = [accuracies, accuracy];
        ps = [ps, stats.p];
    end

    rois(mask_idx).filename = filenames{mask_idx};
    rois(mask_idx).name = masknames{mask_idx};
    rois(mask_idx).accuracies = accuracies;
    rois(mask_idx).ps = ps;

    % some hypothesis testing
    %

    % H0: none of the subjects are significant = all of the correct classificiations were by chance
    %
    n = size(inputs, 1) * numel(subjs);
    k = mean(accuracies) / 100 * n;
    p = 1 - binocdf(k, n, 1/3);
    rois(mask_idx).bino.n = n; 
    rois(mask_idx).bino.k = k; 
    rois(mask_idx).bino.p = p; 

    % H0: none of the subjects are significant = all accuracies should be at 33%  (Kriegeskorte & Bandettini 2007)
    %
    [h,p,ci,stats] = ttest(accuracies/100, 1/3);
    rois(mask_idx).ttest.stats = stats; 
    rois(mask_idx).ttest.ci = ci; 
    rois(mask_idx).ttest.p = p; 
end

filename = sprintf('classify_contrast_%d_%s_%s_%s_%s.mat', glmodel, contrast, method, what, event);
fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), 'rois', 'targets', 'which_rows', 'p', 'alpha', 'Dis', 'Num', 'r', 'direct', 'method', 'runs', 'trials', 'subjs', 'predict_what', 'z_score', 'use_tmaps', 'use_nosmooth');
