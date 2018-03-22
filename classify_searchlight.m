function [Searchlight, idx, x, y, z] = classify_searchlight(start_idx, end_idx, r, event)

% Run classifier for different searchlight spheres
% Similar to rdms_searchlight.m and rdms_get_searchlight.m
%
% INPUT:
% start_idx, end_idx = range of voxel indices where to center the spheres
%                      (those are randomly shuffled so it's just for batching)
% r = radius of sphere, in the coordinate space of the group-level mask
% idx = optionally, list of voxel indices = voxel order. If not supplied, voxels are randomly shuffled.
%
% OUTPUT:
% table_Rho = searchlights x models matrix of correlation coefficients
% table_P = searchlights x models matrix of p-values
% all_subject_rhos = searchlights x models x n_subjects matrix of
%                    correlation coefficients for each subject
% idx = indices of voxels used as sphere centers
% x, y, z = sphere centers in coordinate space of group-level mask

dirname = 'classify';

%method = 'cvglmnet';
%method = 'patternnet';
method = 'cvpatternnet';
runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none'; 

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only


%% Get the searchlight betas
%
whole_brain = load_mask('masks/mask.nii');
[x, y, z] = ind2sub(size(whole_brain), find(whole_brain)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space

% get voxel order based on ANOVA
%
idx = classify_anova_get_voxel_idx(event);
assert(numel(idx) == numel(x));

% reorder voxels
x = x(idx);
y = y(idx);
z = z(idx);

% take the range specified by the user
end_idx = min(end_idx, numel(x));
i = start_idx:end_idx;
x = x(i);
y = y(i);
z = z(i);
idx = idx(i);
disp(end_idx);

% get betas
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

whole_brain_activations = get_activations('masks/mask.nii', event, data, metadata, use_nosmooth);

% get coordinate limits
%
[mask, Vmask] = load_mask('masks/mask.nii');
Vmask.fname = 'masks/searchlight.nii'; % !!!! IMPORTANT! in case we save it

[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);

disp('Classifying...');
tic

% Actually run the classifier
%
for i = 1:numel(x)
    fprintf('Row %d: vox_idx = %d [%d %d %d]\n', i, idx(i), x(i), y(i), z(i));

    % Build spherical mask
    %
    sphere_mask = create_spherical_mask_helper(mask, x(i), y(i), z(i), r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);

    % Get mask betas
    %
    sphere_activations = get_activations_submask(sphere_mask, whole_brain_activations);

    % if a voxel is NaN even for 1 trial, ignore it
    %
    good_voxels = sum(isnan(sphere_activations(which_rows,:)), 1) == 0; 
    if sum(good_voxels) == 0
        assert(use_nosmooth); % doesn't happen if we have smoothing b/c we are using the mask that by definition contains no NaNs
        warning(sprintf('Skipping sphere at %d [%d %d %d] -- no good voxels', idx(i), x(i), y(i), z(i)));
        continue;
    end
    if sum(good_voxels) < size(sphere_activations, 2)
        assert(use_nosmooth);
        sphere_activations = sphere_activations(:,good_voxels);
        warning(sprintf('Sphere at %d [%d %d %d] has only %d good voxels', idx(i), x(i), y(i), z(i), sum(good_voxels)));
    end

    % one classifier for all subjects -- not standard
    %
    %[inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs, sphere_activations, predict_what, z_score, data, metadata);
    %[classifier, outputs, accuracy] = classify_train_helper(method, inputs, targets, runs, trials, subjs, []);

    % one classifier per subject
    %
    accuracies = [];
    ps = [];
    for subj = subjs
        fprintf('    subj %d\n', subj);

        [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, [subj], sphere_activations, predict_what, z_score, data, metadata);
        [classifier, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, [subj], []);

        accuracies = [accuracies, accuracy];
        ps = [ps, stats.p];
    end

    % H0: none of the subjects are significant = all accuracies should be at 33%  (Kriegeskorte & Bandettini 2007)
    %
    [h,p,ci,stats] = ttest(accuracies/100, 1/3);

    % Save results
    %
    mni = cor2mni([x(i) y(i) z(i)], Vmask.mat); % coords in MNI space
    %Searchlight(i).accuracy = accuracy;
    Searchlight(i).accuracies = accuracies;
    Searchlight(i).ps = ps;
    Searchlight(i).p = p;
    Searchlight(i).tstat = stats.tstat;
    Searchlight(i).name = ['sphere_', sprintf('%d_%d_%d', mni), '_', event(1)];
    Searchlight(i).center = [x(i) y(i) z(i)];
    Searchlight(i).radius = r;
    Searchlight(i).center_mni = mni;
end

disp('Classified.');
toc

%% Save output
%

filename = sprintf('searchlight_classifier_%s_%d-%d-%s.mat', method, start_idx, end_idx, event);
fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), 'Searchlight', 'event', 'x', 'y', 'z', 'r', 'idx', 'targets', 'which_rows', 'method', 'runs', 'trials', 'subjs', 'predict_what', 'z_score', 'use_tmaps', 'use_nosmooth');
