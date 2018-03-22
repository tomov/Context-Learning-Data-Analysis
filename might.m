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

disp('getting inputs');
tic
[inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs(1), betas, predict_what, z_score, data, metadata);
toc


[~, labels] = max(targets, [], 2); % from one-hot vector to indices
labelsGroup = data.runId(which_rows); % folds = runs (assumes 1 subject, and that GNB takes prior into account)

classifier = 'gnb_searchmight'; % fast GNB

disp('running searchmight');
tic 
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours);
toc

