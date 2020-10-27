%function [table_Rho, table_P, all_subject_rhos, idx, x, y, z] = rdms_searchlight(mni, r)
% copy of rdms_searchlight_isl.m

UNFINISHED


dirname = 'rdms';
filename = sprintf('rdms_searchlight_isl.mat');
r = 4; % mm

masks = isl_create_masks(false, r);

use_tmaps = false;
use_nosmooth = false;

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only


[Model, control_model_idxs, params] = rdms_get_model_isl(data, metadata, which_rows);

fprintf('SAVING %s\n', filename);
save(fullfile(dirname, filename), '-v7.3');

%load(fullfile(dirname, filename));


ROIs = rdms_get_rois_isl(data, metadata, which_rows, x, y, z, r, true, use_tmaps, use_nosmooth); % use pregen'd betas, use tmaps, use nosmooth


%% Get second-order RDM
%
assert(isequal(Model(control_model_idxs(1)).name, 'time'));
assert(isequal(Model(control_model_idxs(2)).name, 'run'));
rows = ROIs;
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
names = {ROIs.name};
Ts = [];
Ps = [];
for r = 1:length(ROIs)
    rho1 = squeeze(all_subject_rhos(r,1,:));
    for m = 2:4
        rho2 = squeeze(all_subject_rhos(r,m,:));
        [h,p,ci,stat] = ttest(atanh(rho1), atanh(rho2));
        Ts(r,m) = stat.tstat;
        Ps(r,m) = p;
    end
end

tbl = table(names', Ts, Ps)
