%function [table_Rho, table_P, all_subject_rhos, idx, x, y, z] = rdms_searchlight(mni, r)
% copy of rdms_searchlight.m

% TODO dedupe w/ isl_create_masks
mni = [ 44 20 -6; ... % AI  context
        -46 -38 44; ... % IPL  context
        42 24 24; ... % IFG  Agency
        36 16 6; ... % AI  agency
        42 28 26; ... % IFG  vgdl
        -30 28 2; ... % AI vgdl
        -50 44 12; ... % IFG  vgdl
        ];

r = 4 / 1.5; % mm -> voxels

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

dirname = 'rdms';
filename = sprintf('rdms_searchlight_isl.mat');

%{
use_tmaps = false;
use_nosmooth = false;

[mask, Vmask] = load_mask('masks/mask.nii');
cor = mni2cor(mni, Vmask.mat);
x = cor(:,1);
y = cor(:,2);
z = cor(:,3);

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only


[Model, control_model_idxs, params] = rdms_get_model_isl(data, metadata, which_rows);

fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), '-v7.3');
%}
load(fullfile(dirname, filename));


Searchlight = rdms_get_searchlight(data, metadata, which_rows, x, y, z, r, true, use_tmaps, use_nosmooth); % use pregen'd betas, use tmaps, use nosmooth


%% Get second-order RDM
%
assert(isequal(Model(control_model_idxs(1)).name, 'time'));
assert(isequal(Model(control_model_idxs(2)).name, 'run'));
rows = Searchlight;
cols = Model;

[table_Rho, table_H, table_T, table_P, all_subject_rhos] = rdms_second_order(metadata, rows, cols, control_model_idxs, false, [], []);

%% Save output
%


fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), '-v7.3');
%save(fullfile(dirname, filename), 'Model', 'table_Rho', 'table_T', 'table_P', 'all_subject_rhos', 'x', 'y', 'z', 'events', 'mni', 'r', 'idx', 'params', 'use_tmaps', 'use_nosmooth', 'which_rows');



% paired t-tests
% MCMC vs. ideal observer
%
names = {Searchlight.name};
Ts = [];
Ps = [];
for r = 1:length(Searchlight)
    rho1 = squeeze(all_subject_rhos(r,1,:));
    for m = 2:4
        rho2 = squeeze(all_subject_rhos(r,m,:));
        [h,p,ci,stat] = ttest(atanh(rho1), atanh(rho2));
        Ts(r,m) = stat.tstat;
        Ps(r,m) = p;
    end
end


ms = {Model.name};
ms(1:4)
table(names', table_T(:,1:4), table_P(:,1:4), 'VariableNames', {'ROI', 'RSA t-stat (vs. 0)', 'RSA p-value'})


ms(2:4)
tbl = table(names', Ts(:,2:4), Ps(:,2:4), 'VariableNames', {'ROI', 'paired t-stat (vs. ideal)', 'p-value'})
