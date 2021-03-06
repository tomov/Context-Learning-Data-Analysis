% Display Fmap computed by classify_anova.m
%

%event = 'feedback_onset';
event = 'trial_onset';

load(['results/classify_anova_', event, '.mat']);

[mask, V] = load_mask('masks/mask.nii');
V.fname = fullfile('results', ['classify_anova_', event, '_Fmap.nii']);
[x, y, z] = ind2sub(size(mask), find(mask));

assert(numel(x) == numel(F_values));

Fmap = nan(size(mask));
Fmap(mask) = F_values;

spm_write_vol(V, Fmap);

bspmview(V.fname, 'masks/mean.nii');
