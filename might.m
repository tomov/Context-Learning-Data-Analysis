%{
addpath('/Users/momchil/Dropbox/Research/libs/SearchmightToolbox.Darwin_i386.0.2.5');

maskfile = 'masks/mask.nii';
event = 'feedback_onset';
use_nosmooth = true;
r = 2.6667;

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

% fix neighbors according to spherical searchlight
%
[meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(mask, meta.colToCoord, r);

betas = get_betas(maskfile, event, data, metadata, use_nosmooth);
dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
nvoxels = size(betas, 2);

classifier = 'gnb_searchmight'; % fast GNB

which_rows = data.which_rows & data.isTrain & strcmp(data.participant, 'con001');

%}

%{
disp('getting inputs');
tic
[inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs(1), betas, predict_what, z_score, data, metadata);
toc


[~, labels] = max(targets, [], 2); % from one-hot vector to indices
% TODO compare w/ balanced folds
labelsGroup = data.runId(which_rows); % folds = runs (assumes 1 subject, and that GNB takes prior into account)

classifier = 'gnb_searchmight'; % fast GNB

disp('running searchmight');
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours);
%}

dirname = 'might';

filename = sprintf('%s_accuracy_subj=%d_folds=%d_r=%.4f_use_nosmooth=%d.nii', classifier, subjs(1), max(labelsGroup), r, use_nosmooth);

% initialize an empty map
[~, V, amap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
amap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write map
amap(mask) = am * 100;
spm_write_vol(V, amap);

% visualize
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);
