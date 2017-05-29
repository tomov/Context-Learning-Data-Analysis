function betas = load_betas_trial_onset(mask, data, metadata)

% Load all the betas corresponding to activations at trial onset for the
% given mask. Does this for all trials for all runs for all subjects.
%
% INPUT:
% mask = path to .nii file with the mask, e.g. 'masks/hippocampus.nii'
% data, metadata = subject behavioral data as output by load_data
%
% OUTPUT:
% betas = [nRows x nVoxels] beta coefficients for the mask for each trial.
%         Each row of this matrix corresponds to a single trial == a single 
%         row in data.
%

EXPT = context_expt();
glmodel = 60; % this is the one that has a regressor for each trial onset
subjs = getGoodSubjects();
runs = 1:metadata.runsPerSubject;
trials = 1:metadata.trialsPerRun;

betas = [];

for subj = subjs
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    
    for run = runs
        for trial = trials
            regressor = ['Sn(', num2str(run), ') trial_onset_', num2str(trial)];
            
            beta = ccnl_get_beta(EXPT, glmodel, regressor, mask, subj);
            assert(size(beta, 1) == 1);
            %beta(isnan(beta)) = 0; 
            
            which_rows = strcmp(data.participant, metadata.allSubjects{subj}) & data.runId == run & data.newTrialId == trial;
            assert(sum(which_rows) == 1);
            
            % initialize betas array if we haven't yet
            %
            if isempty(betas)
                betas = nan(size(data.participant, 1), size(beta, 2));
            end
            % append beta
            betas(which_rows, :) = beta;
            
            fprintf('loaded subj %d, run %d, trail %d (%s): %d voxels\n', subj, run, trial, regressor, sum(~isnan(beta)));
        end
    end
end
