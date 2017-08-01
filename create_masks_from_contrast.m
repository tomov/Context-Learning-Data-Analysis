function mask_filenames = create_masks_from_contrast(EXPT, model, contrast, p)
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction and create masks from them. Uses the same logic as bspmview,
% and as a sanity check prints the results table -- should be the same as
% the one from bspmview.
%
% INPUT:
% EXPT = experiment structure, e.g. context_expt()
% model = GLM number, e.g. 154
% contrast = contrast, e.g. 'KL_weights - KL_structures'
% p = optional p-value threshold; defaults to 0.001
%
% OUTPUT:
% mask_filenames = the masks
%

if ~exist('p', 'var')
    p = 0.001;
end

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
[C, CI, region, extent, stat, mni, cor, results_table] = bspm_extract_clusters(spmT, p, '+/-');
[~, V, Y] = load_mask(spmT);
V.fname = 'temp/temp.nii'; % <-- change immediately so we don't overwrite it accidentally

disp(results_table);

mask_filenames = [];

% create a mask with each cluster
%
for i = 1:size(region, 1)
    
    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);
    
    clust_idx = CI(x,y,z);
    mask = CI == clust_idx;
    
    filename = sprintf('glm%d %s cluster t=%.3f extent=%d roi=%s peak=[%d %d %d].nii', model, contrast, stat(i), extent(i), region{i}, mni(i,1), mni(i,2), mni(i,3));
    V.fname = fullfile('masks', filename);
    disp(V.fname);
    spm_write_vol(V, mask);
    
    mask_filenames{i} = V.fname;
end