% For each voxel, see if it varies significantly with condition (ANOVA), 
% get F-statistic and p-value to obtain a F-map
%

use_tmaps = false;
use_nosmooth = false;
if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end
event = 'feedback_onset';

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% get betas
%
%betas = get_activations('masks/mask.nii', event, data, metadata, use_nosmooth);

which_rows = data.which_rows & data.isTrain;

condition_labels = containers.Map(metadata.conditions, {1, 2, 3});
conds = values(condition_labels, data.contextRole);
conds = cell2mat(conds);


for vox_idx = 1:size(betas, 2) 
    [p, anovatab, stats] = anova1(betas(which_rows, vox_idx), conds(which_rows), 'off');
    F = anovatab{2,5};

    p_values(vox_idx) = p;
    F_values(vox_idx) = F;

    fprintf('%d: F = %.5f, p = %.5f\n', vox_idx, F, p);
end

outfilename = fullfile('results', sprintf('classify_anova_%s.mat', event));

save(outfilename, 'p_values', 'F_values', 'use_tmaps', 'use_nosmooth', 'event', 'which_rows');
