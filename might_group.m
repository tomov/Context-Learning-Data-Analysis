% Redo group-level analysis of searchmight results & visualize them
%


dirname = 'might';

maskfile = 'masks/mask.nii';
event = 'trial_onset';
r = 2.6667;

use_tmaps = false;
use_nosmooth = true; 

runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none'; 

classifier = 'lda_shrinkage';
%classifier = 'gnb_searchmight';


alpha = 0.05; % this parameter is used in redoing the stats here

[mask] = load_mask(maskfile);

accuracy_map_template = '%s_accuracy_%s_subj=%d_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii';
pvalue_map_template = '%s_p-value_%s_subj=%d_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii';

nvoxels = sum(mask(:));
howmany = zeros(nvoxels, 1); % for each voxel, how many subjects have it significant 

for subj = subjs 
    fprintf('    subj %d\n', subj);

    filename = sprintf(accuracy_map_template, classifier, event, subj, 3, r, z_score, use_nosmooth, use_tmaps);
    [~,~,amap] = load_mask(fullfile(dirname, filename));

    filename = sprintf(pvalue_map_template, classifier, event, subj, 3, r, z_score, use_nosmooth, use_tmaps);
    [~,~,pmap] = load_mask(fullfile(dirname, filename));

    am = amap(mask);
    pm = pmap(mask);

    % Benjamini-Hochberg FDR
    %
    [fdr] = mafdr(pm, 'BHFDR', true);
    howmany = howmany + (fdr < alpha);

    % Storey (2002) pFDR
    %
    %[pfdr, qm] = mafdr(pm);
    %howmany = howmany + (qm < alpha);
end

[~, V, countmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
V.fname = 'temp.nii'; % change immediately!
countmap(:) = NaN;
countmap(mask) = howmany;
spm_write_vol(V, countmap);

bspmview(V.fname, 'masks/mean.nii');
