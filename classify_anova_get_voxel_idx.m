function [idx] = classify_anova_get_voxel_idx(filename)
% Helper function to get the voxel indices in descending order of F-values (after ANOVA looking for condition-sensitive voxels)
% i.e. voxels most likely to correspond to condition come first
%

if ~exist('filename', 'var') || isempty(filename)
    % by default, load voxels sensitive to condition @ feedback onset
    event = 'feedback_onset';
    filename = fullfile('results', ['classify_anova_', event, '.mat']);
end

load(filename);

[~, idx] = sort(F_values, 'descend');
