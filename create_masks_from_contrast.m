function [mask_filenames, mask_names] = create_masks_from_contrast(EXPT, model, contrast, p, direct, alpha, Dis, Num)
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction and create masks from them. Uses the same logic as bspmview,
% and as a sanity check prints the results table -- should be the same as
% the one from bspmview.
%
% INPUT:
% EXPT = experiment structure, e.g. context_expt()
% model = GLM number, e.g. 154
% contrast = contrast, e.g. 'KL_weights - KL_structures'
% p = p-value threshold; defaults to 0.001
% direct = sign of the activations; should be one of +, -, or +/-;
%          defaults to +/-
% alpha = significance level for cluster FWE correction, default 0.05 in bspmview
% Dis = separation for cluster maxima, default 20 in bspmview
% Num = numpeaks for cluster maxima, default 3 in bspmview
%
% OUTPUT:
% mask_filenames = paths to the '.nii' files; saved by default in 'masks'
%                  subdirectory
% mask_names = shorter names of the masks to be used e.g. in RDMs
%

assert(ismember(direct, {'+/-', '+', '-'}));

% extract the clusters
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, model, contrast, p, direct, alpha, Dis, Num);

mask_filenames = [];
mask_names = [];

% create a mask with each cluster
%
for i = 1:size(region, 1)
    
    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);
    
    clust_idx = CI(x,y,z);
    mask = CI == clust_idx;
    
    filename = sprintf('glm%d_%s_cluster_t=%.3f_extent=%d_roi=%s_peak=[%d_%d_%d].nii', model, contrast, stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));
    V.fname = fullfile('masks', filename);
    disp(V.fname);
    spm_write_vol(V, mask);
    mask_filenames{i} = V.fname;
    
    maskname = sprintf('GLM%d_%s_%s_%.3f_%d', model, contrast, region{i}, stat(i), extent(i));
    mask_names{i} = maskname;
end
