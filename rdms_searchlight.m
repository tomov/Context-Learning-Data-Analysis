function [table_Rho, table_P, all_subject_rhos, idx, x, y, z] = rdms_searchlight(start_idx, end_idx, r)

% Compare RDMs from different searchlight spheres with model RDMs
% Similar to rdms.m
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

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only

%% Get the searchlight RDMs
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
idx = start_idx:end_idx;
x = x(idx);
y = y(idx);
z = z(idx);

Searchlight = rdms_get_searchlight(data, metadata, which_rows, x, y, z, r, true, true, false); % use pregen'd betas, use tmaps, use nosmooth

%% Get the model RDMs
%
Model = rdms_get_model(data, metadata, which_rows);


%% Get second-order RDM
%
control_model_idxs = [8, 12]; % #KNOB control for time and run
assert(isequal(Model(8).name, 'time'));
assert(isequal(Model(12).name, 'run'));
rows = Searchlight;
cols = Model;

[table_Rho, table_H, table_T, table_P, all_subject_rhos] = rdms_second_order(metadata, rows, cols, control_model_idxs, false, [], []);


%% Save output
%

%which = table_Rho > 0 & table_P < 0.05 / numel(table_P); % Bonferroni correction

filename = sprintf('searchlight_%d-%d.mat', start_idx, end_idx);
save(fullfile('rdms', filename), 'table_Rho', 'table_T', 'table_P', 'all_subject_rhos', 'x', 'y', 'z', 'r', 'idx');
