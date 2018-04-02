function run_betas = load_run_betas(glmodel, regressor_prefix, voxels, is_mni)
% For a given GLM, get the single-subject betas for each run for a given
% regressor and a set of voxels.
% This is used e.g. to get the KL divergence betas so then you can do a
% within-subject correlation with behavior on the test trials
%
% INPUT:
% glmodel = GLM from which to take the betas
% regressor_prefix = prefix of the regressor of interest
% voxels = V x 3 list of [x y z] voxel coordinates in MNI space (V = # of
%          voxels)
% is_mni = optional whether voxels are in MNI or native space (default true)
%
% OUTPUT:
% run_betas = N x R x V vector of betas, where N = # subjects, R = # runs
%             per subject, and V = # of voxels)
%

EXPT = context_expt();

if ~exist('is_mni', 'var')
    is_mni = true;
end

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

        % For each run, get the betas for all voxels
        %
        beta_idx = regressors_idxs(run);
        assert(~isempty(strfind(SPM.xX.name{beta_idx}, regressor_prefix)));
        V = spm_vol(fullfile(subjdir, sprintf('beta_%04d.nii', beta_idx)));
        Y = spm_read_vols(V);

        % optionally convert from MNI to native coordinates
        if is_mni
            cors = mni2cor(voxels, V.mat);
        else
            cors = voxels;
        end 
        vox_idx = sub2ind(size(Y), cors(:,1), cors(:,2), cors(:,3)); % get voxels as indices for conveniencet

        run_betas(subj_idx, run, :) = Y(vox_idx);

        % now get the actual voxels values <-- lame
        %for v_idx = 1:size(voxels, 1)
        %    cor = cors(v_idx, :);
        %    value = Y(cor(1), cor(2), cor(3));
        %    run_betas(subj_idx, run, v_idx) = value;
        %end
    end
end

