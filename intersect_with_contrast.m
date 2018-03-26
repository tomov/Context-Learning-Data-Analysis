function intersect_with_contrast(infile, outfile, EXPT, model, contrast, p, direct, alpha, Dis, Num)
% Intersect a .nii file with a GLM contrast.
%
% INPUT:
% infile = path to .nii file 
% outfile = path to where to save the output .nii file
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
%


assert(ismember(direct, {'+/-', '+', '-'}));

[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, model, contrast, p, direct, alpha, Dis, Num);

[~, Vin, Yin] = load_mask(infile);

Vout = Vin; 
Vout.fname = outfile; % change immediately!
Yout = nan(size(Yin));

for i = 1:size(region, 1)    

    clust_idx = CI(cor(i,1), cor(i,2), cor(i,3));
    Yout(CI == clust_idx) = Yin(CI == clust_idx);
end

spm_write_vol(Vout, Yout);
bspmview(Vout.fname, 'masks/mean.nii');
