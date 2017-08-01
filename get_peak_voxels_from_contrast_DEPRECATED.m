function [mni, region] = get_peak_voxels_from_contrast(EXPT, model, contrast, p, direct)
% TODO dedupe with extract_clusters
%
% Given a contrast, extract all the peak voxels from the clusters after cluster FWE
% correction and create masks from them. Uses the same logic as bspmview,
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
% mni = MNI coordinates of peak voxels
% region = labels for the peak voxels

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


