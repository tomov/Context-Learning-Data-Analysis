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

idx = 1:100;

Searchlight = rdms_get_searchlight(data, metadata, which_rows, x(idx), y(idx), z(idx), r);
showRDMs(Searchlight, 1);

%% Get the model RDMs
%
Model = rdms_get_model(data, metadata, which_rows);
showRDMs(Model, 2);

...

    rows = Neural;
    cols = Model;
    [table_Rho, table_H, table_P, all_subject_rhos, lmes] = rdms_second_order(metadata, rows, cols, control_model_idxs, do_LME, lme_neural_idxs, lme_model_idxs);


