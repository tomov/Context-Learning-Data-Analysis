function [rois, targets, which_rows] = classify_contrast(EXPT, glmodel, contrast, what)

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

if ~exist('contrast', 'var')
    contrast = regressor;
end

use_tmaps = false;
use_nosmooth = false;

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
direct = '+';

% classifier params
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
%method = 'cvglmnet';
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
event = 'feedback_onset';

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

    [m, V] = load_mask(filenames{mask_idx});
    activations = get_activations_submask(m, whole_brain_activations);

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

    % H0: none of the subjects are significant
    n = size(inputs, 1) * numel(subjs);
    k = mean(accuracies) / 100 * n;
    p = 1 - binocdf(k, n, 1/3);

    save shit.mat;

    rois(mask_idx).filename = filenames{mask_idx};
    rois(mask_idx).name = masknames{mask_idx};
    rois(mask_idx).accuracies = accuracies;
    rois(mask_idx).ps = ps;
    rois(mask_idx).p = p; 

    fprintf(' ROI p = %f\n', p);
end

filename = sprintf('classify_contrast_%d_%s_%s.mat', glmodel, contrast, what);
fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), 'rois', 'targets', 'which_rows', 'p', 'alpha', 'Dis', 'Num', 'r', 'direct', 'method', 'runs', 'trials', 'subjs', 'predict_what', 'z_score');
