function [mask_filename, mask_name] = create_mask_from_contrast(EXPT, model, contrast, p, direct, alpha, Dis, Num)
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction and create a single mask from them. Essentially the same ask create_masks_from_contrast,
% except it creates a single mask rather than a mask for each cluster
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
% mask_filename = path to the '.nii' file; saved by default in 'masks'
%                  subdirectory
% mask_name = shorter name of the mask to be used e.g. in RDMs
%

assert(ismember(direct, {'+/-', '+', '-'}));

% extract the clusters
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, model, contrast, p, direct, alpha, Dis, Num);

mask = logical(zeros(size(Y)));

% create a mask with each cluster
%
for i = 1:size(region, 1)
    
    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);
    
    clust_idx = CI(x,y,z);
    mask = mask | (CI == clust_idx);
end

filename = sprintf('glm%d_%s.nii', model, contrast);
V.fname = fullfile('masks', filename);
disp(V.fname);
spm_write_vol(V, mask);
mask_filename = V.fname;

maskname = sprintf('GLM%d_%s', model, contrast);
mask_name = maskname;
