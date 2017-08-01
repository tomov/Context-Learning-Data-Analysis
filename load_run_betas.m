function run_betas = load_run_betas(glmodel, regressor_prefix, voxels)
% For a given GLM, get the single-subject betas for each run for a given
% regressor.
% This is used e.g. to get the KL divergence betas so then you can do a
% within-subject correlation with behavior on the test trials
%
% INPUT:
% glmodel = GLM from which to take the betas
% regressor_prefix = prefix of the regressor of interest
% voxels = V x 3 list of [x y z] voxel coordinates in MNI space (V = # of
%          voxels)
%
% OUTPUT:
% run_betas = N x R x V vector of betas, where N = # subjects, R = # runs
%             per subject, and V = # of voxels)
%

EXPT = context_expt();

% Load behavior
%
[~, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();


run_betas = nan(metadata.N, metadata.runsPerSubject, size(voxels, 1));

% Iterate over subjects
%
subj_idx = 0; % 1..20
for subj = goodSubjects % 1..25
    subj_idx = subj_idx + 1;
    fprintf('subj %d (idx %d)\n', subj, subj_idx);

    % get the KL regressor idxs, so we can get the betas from our peak voxels 
    %
    subjdir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(subjdir,'SPM.mat'));
    regressors_idxs = find(~cellfun('isempty', strfind(SPM.xX.name, regressor_prefix)));
    assert(numel(regressors_idxs) == metadata.runsPerSubject);

    % Iterate over runs
    %
    for run = 1:metadata.runsPerSubject

        % For each run, get the beta for each ROI peak voxel
        %
        beta_idx = regressors_idxs(run);
        assert(~isempty(strfind(SPM.xX.name{beta_idx}, regressor_prefix)));
        % this does not give you the full 3D volume, just the non-nan
        % betas in a linear vector.
        V = spm_vol(fullfile(subjdir, sprintf('beta_%04d.nii', beta_idx)));
        Y = spm_read_vols(V);
        % now get the actual voxels values
        for v_idx = 1:size(voxels, 1)
            voxel = voxels(v_idx, :);
            cor = mni2cor(voxel, V.mat);
            value = Y(cor(1), cor(2), cor(3));
            run_betas(subj_idx, run, v_idx) = value;
        end
    end
end

