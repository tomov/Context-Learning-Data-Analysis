function [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, model, contrast, p, direct)
% Wrapper around bspm_extract_clusters
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction. Uses the same logic as bspmview,
% and as a sanity check prints the results table -- should be the same as
% the one from bspmview.
%
% INPUT:
% EXPT = experiment structure, e.g. context_expt()
% model = GLM number, e.g. 154
% contrast = contrast, e.g. 'KL_weights - KL_structures'
% p = optional p-value threshold; defaults to 0.001
% direct = optional sign of the activations; should be one of +, -, or +/-;
%          defaults to +/-
%
% OUTPUT:
% V = SPM volume of the t-map, with the filename changed so we don't
%     overwrite it accidentally
% Y = the actual t-map
% C = volume with cluster size for each voxel
% CI = volume with cluster index for each voxel <-- that's the name of the
%      game; 
% region = labels for the peak voxels
% extent = size of the cluster
% stat = t-statistic for the peak voxels
% mni = MNI coordinates of peak voxels
% cor = coordinates of peak voxel in native space (can be used as indices
%       in the C and CI volumes)
% results_table = what Save Results Table in bspmview would spit out 
% 


if ~exist('p', 'var')
    p = 0.001;
end
if ~exist('direct', 'var')
    direct = '+/-';
end
assert(ismember(direct, {'+/-', '+', '-'}));

% find the contrast
%
modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
load(fullfile(modeldir,'contrasts'));
ix = find(strcmp(contrasts,contrast));
if isempty(ix)
    error('Contrast not found');
end
spmT = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(ix)],'spmT_0001.nii');

% extract the clusters
%
[C, CI, region, extent, stat, mni, cor, results_table] = bspm_extract_clusters(spmT, p, direct);
[~, V, Y] = load_mask(spmT);
V.fname = 'temp/temp.nii'; % <-- change immediately so we don't overwrite it accidentally