addpath('/Users/momchil/Dropbox/Research/libs/SearchmightToolbox.Darwin_i386.0.2.5');

dirname = 'might';

maskfile = 'masks/mask.nii';
%event = 'feedback_onset';
event = 'trial_onset';
r = 2.6667;

use_tmaps = true;
use_nosmooth = true;

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end

runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
%z_score = 'z-run'; 
z_score = 'z-none'; 


[mask] = load_mask(maskfile);
% load instead
%
%[meta] = createMetaFromMask(mask); <--- takes forever; load it instead (precomputed)
load('might/meta.mat');

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

activations = get_activations(maskfile, event, data, metadata, use_nosmooth);
dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
nvoxels = size(activations, 2);

if use_tmaps
    % some NaN's become 0's in the t-maps
    % note that the betas (and the t-values) can never really be 0
    %
    activations(activations == 0) = NaN;
end

% fix neighbors according to spherical searchlight
%
[meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(meta.colToCoord, r, mask, activations, data.which_rows & data.isTrain);





ams = nan(numel(subjs), nvoxels); % accuracy map for each subject 

alpha = 0.05; % for pFDR q-values
howmany = zeros(nvoxels, 1); % for each voxel, how many subjects have it significant (according to pFDR)

subj_idx = 0;
for subj = subjs 
    subj_idx = subj_idx + 1;
    fprintf('    subj %d\n', subj);

    [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subj, activations, predict_what, z_score, data, metadata);

    [~, labels] = max(targets, [], 2); % from one-hot vector to indices
    % TODO compare w/ balanced folds
    %labelsGroup = data.runId(which_rows); % folds = runs; WRONG -- GNB does not take prior properly into account; need balanced folds
    [labelsGroup, kfolds] = balanced_folds(runs, subj, trials, targets); % 3 balanced folds
    assert(kfolds == 3);

    classifier = 'gnb_searchmight'; % fast GNB; quick-n-dirty; gives weird axial dashes in accuracy maps, probably b/c of the way the scanner takes the images
    %classifier = 'lda_shrinkage'; % based on Pereira and Botvinick (2010)

    disp('running searchmight');
    [am,pm] = computeInformationMap(inputs,labels,labelsGroup,classifier,'searchlight', ...
                                    meta.voxelsToNeighbours,meta.numberOfNeighbours);
    ams(subj_idx,:) = am;

    [~, qm] = mafdr(pm'); % Storey (2002)
    howmany = howmany + (qm < alpha);

    disp('saving ouput');
    filename = sprintf('%s_accuracy_%s_subj=%d_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, subj, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
    % initialize an empty accuracy map
    [~, V, amap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
    amap(:) = NaN; % clear
    V.fname = fullfile(dirname, filename); % change immediately!

    % write accuracy map
    amap(mask) = am * 100;
    spm_write_vol(V, amap);

    % visualize
    %struc = fullfile('masks','mean.nii');
    %bspmview(V.fname, struc);
end


%% write map w/ # subjects for which voxel is significant (with pFDR) based on Pereira & Botvinick 2010
%

% initialize empty countmap 
filename = sprintf('%s_accuracy_countmap_%s_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
[~, V, countmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
countmap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write countmap
%
[h, p, ci, stats] = ttest(ams, 1/3);
countmap(mask) = howmany;
spm_write_vol(V, countmap);

% visualize countmap
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);




%% write t-map based on Kriegeskorte & Bandettini 2007 
% note those have lots of negative t-values b/c "chance" is not 1/3 according to the classifier but slightly below it
% so don't use em
%

%{
% initialize empty tmap
filename = sprintf('%s_accuracy_tmap_%s_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
[~, V, tmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
tmap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write tmap
% for each voxel, t-test subject accuracies against chance 
%
[h, p, ci, stats] = ttest(ams, 1/3);
tmap(mask) = stats.tstat;
spm_write_vol(V, tmap);

% visualize tmap
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);

%}
