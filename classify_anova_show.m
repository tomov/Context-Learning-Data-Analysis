% Display Fmap computed by classify_anova.m
%

event = 'feedback_onset';

load(['results/classify_anova_', event, '.mat']);

[mask, V] = load_mask('masks/mask.nii');
V.fname = fullfile('results', ['classify_anova_', event, '_Fmap.nii']);
[x, y, z] = ind2sub(size(mask), find(mask));

assert(numel(x) == numel(F_values));

Fmap = nan(size(mask));
Fmap(mask) = F_values;
%for i = 1:numel(x)
%    Fmap(x(i), y(i), z(i)) = F_values(i);
%end

spm_write_vol(V, Fmap);

bspmview(V.fname, 'masks/mean.nii');
