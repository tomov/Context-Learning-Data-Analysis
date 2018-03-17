function [idx] = classify_anova_get_voxel_idx(event)
% Helper function to get the voxel indices in descending order of F-values (after ANOVA looking for condition-sensitive voxels)
% i.e. voxels most likely to correspond to condition come first
%

assert(ismember(event, {'feedback_onset', 'trial_onset'}));

filename = fullfile('results', ['classify_anova_', event, '.mat']);
load(filename);

[~, idx] = sort(F_values, 'descend');
