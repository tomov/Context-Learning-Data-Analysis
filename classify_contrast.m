function [rois, targets, which_rows] = classify_contrast(glmodel, contrast, what)

% Train classifier using ROIs from a contrast.
%
% INPUT:
% glmodel = glm as in context_create_multi.m
% contrast = contrast from which to extract the clusters; by
% what = 'sphere', or 'cluster' -- what area to take around the
%        peak voxel from each cluster
%
% OUTPUT:
%

if ~exist('contrast', 'var')
    contrast = regressor;
end

% fmri cluster params
%
EXPT = context_expt();
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
method = 'cvglmnet';
runs = 1:metadata.runsPerSubject;
trials = 6:metadata.trainingTrialsPerRun; % TODO arbitrary
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'none';
event = 'feedback_onset';


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
    fprintf('Cluster %d: mask %d\n', mask_idx, filenames{mask_idx});

    rois(mask_idx).filename = filenames{mask_idx};
    rois(mask_idx).rdm_name = masknames{mask_idx};
    [classifier, ~, targets, outputs, which_rows, accuracy] = classify_train(method, runs, trials, subjs, filenames{mask_idx}, predict_what, z_score, event);
    rois(mask_idx).classifier = classifier;
    rois(mask_idx).accuracy = accuracy;
    rois(mask_idx).outputs = outputs;
end
