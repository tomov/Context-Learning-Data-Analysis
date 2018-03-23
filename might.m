%{
addpath('/Users/momchil/Dropbox/Research/libs/SearchmightToolbox.Darwin_i386.0.2.5');

dirname = 'might';

maskfile = 'masks/mask.nii';
event = 'feedback_onset';
r = 2.6667;

use_tmaps = false;
use_nosmooth = true;

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_activations;
    load_activations = @load_activations;
end

runs = 1:9; 
trials = 6:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none'; 


[mask] = load_mask(maskfile);
% load instead
%
%[meta] = createMetaFromMask(mask);
load('might/meta.mat');

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

activations = get_activations(maskfile, event, data, metadata, use_nosmooth);
dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
nvoxels = size(activations, 2);

classifier = 'gnb_searchmight'; % fast GNB

% fix neighbors according to spherical searchlight
%
[meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(meta.colToCoord, r, mask, activations, data.which_rows & data.isTrain);

%}


%{
ams = nan(numel(subjs), nvoxels); % accuracy map for each subject

tic
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

    classifier = 'gnb_searchmight'; % fast GNB

    disp('running searchmight');
    [am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                    meta.voxelsToNeighbours,meta.numberOfNeighbours);
    ams(subj_idx,:) = am;

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


%}


% initialize empty tmap
filename = sprintf('%s_accuracy_tmap_%s_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
[~, V, tmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
tmap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write tmap
% for each voxel, t-test subject accuracies against chance (Kriegeskorte & Bandettini 2007)
%
[h, p, ci, stats] = ttest(ams, 1/3);
tmap(mask) = stats.tstat;
spm_write_vol(V, tmap);

% visualize tmap
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);


