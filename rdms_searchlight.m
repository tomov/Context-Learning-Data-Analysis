% Compare RDMs from different searchlight spheres with model RDMs
% Similar to rdms.m
%


%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only

%% Get the searchlight RDMs
%
rng default; % make sure it's the same every time

whole_brain = load_mask('masks/mask.nii');
[x, y, z] = ind2sub(size(whole_brain), find(whole_brain)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space

idx = randperm(length(x));
x = x(idx);
y = y(idx);
z = z(idx);
r = 5;

idx = 201:300;

Searchlight = rdms_get_searchlight(data, metadata, which_rows, x(idx), y(idx), z(idx), r);

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

[table_Rho, table_H, table_P, all_subject_rhos] = rdms_second_order(metadata, rows, cols, control_model_idxs, false, [], []);



which = table_Rho > 0 & table_P < 0.05 / numel(table_P); % Bonferroni correction


%%
%save('rdms/spotlight_101-200.mat', 'table_H', 'table_P', 'table_Rho', 'all_subject_rhos', 'x', 'y', 'z', 'r', 'Model')
