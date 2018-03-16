function [table_Rho, table_P, all_subject_rhos, idx, x, y, z] = classify_searchlight(start_idx, end_idx, r)

% Run classifier for different searchlight spheres
% Similar to rdms_searchlight.m
%
% INPUT:
% start_idx, end_idx = range of voxel indices where to center the spheres
%                      (those are randomly shuffled so it's just for batching)
% r = radius of sphere, in the coordinate space of the group-level mask
%
% OUTPUT:
% table_Rho = searchlights x models matrix of correlation coefficients
% table_P = searchlights x models matrix of p-values
% all_subject_rhos = searchlights x models x n_subjects matrix of
%                    correlation coefficients for each subject
% idx = indices of voxels used as sphere centers
% x, y, z = sphere centers in coordinate space of group-level mask

dirname = 'classify';

method = 'cvglmnet';
runs = 1:9; 
trials = 5:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none'; 

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only


%% Get the searchlight betas
%
rng default; % make sure it's the same every time

whole_brain = load_mask('masks/mask.nii');
[x, y, z] = ind2sub(size(whole_brain), find(whole_brain)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space

% random shuffle them (deterministically)
idx = randperm(length(x));
x = x(idx);
y = y(idx);
z = z(idx);

% take the range specified by the user
end_idx = min(end_idx, numel(x));
idx = start_idx:end_idx;
x = x(idx);
y = y(idx);
z = z(idx);
disp(end_idx);

Searchlight = classify_get_searchlights(data, metadata, which_rows, x, y, z, r, true, false, false); % use pregen'd betas, use tmaps, use nosmooth

disp('Classifying...');

% Actually run the classifier
%
for i = 1:numel(Searchlight)
    fprintf('row %d\n', i);

    [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs, Searchlight(i).activations, predict_what, z_score, data, metadata);

    outFilename = []; % don't save it
    [classifier, outputs, accuracy] = classify_train_helper(method, inputs, targets, runs, trials, subjs, outFilename);

    Searchlight(i).classifier = classifier;
    Searchlight(i).outputs = outputs;
    Searchlight(i).accuracy = accuracy;
end

disp('Classified.');

%% Save output
%

filename = sprintf('searchlight_classifier_%d-%d.mat', start_idx, end_idx);
fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), 'Searchlight', 'x', 'y', 'z', 'r', 'idx', 'params', 'which_structures', 'targets', 'which_rows');
