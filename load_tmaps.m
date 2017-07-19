function tmaps = load_tmaps(mask, regressor_prefix, data, metadata)

% Load all the t-maps corresponding to activations at trial onset or feedback onset for the
% given mask and regressor prefix. Does this for all trials for all runs for all subjects.
% TODO dedupe w/ load_betas
%
% INPUT:
% mask = path to .nii file with the mask, e.g. 'masks/hippocampus.nii'
% regressor_prefix = 'trial_onset' or 'feedback_onset'; assumes regressor
%                    names are of the form {regressor prefix}_{trial id}
%                    and that the contrast names used to generate the
%                    t-maps have the exact same names as the regressor
%                    names (including the Sn(#) prefix)
% data, metadata = subject behavioral data as output by load_data
%
% OUTPUT:
% tmaps = [nRows x nVoxels] t statistics for the voxels in the mask for each trial.
%         Each row of this matrix corresponds to a single trial == a single 
%         row in data.
%

assert(ismember(regressor_prefix, {'trial_onset', 'feedback_onset'}));

EXPT = context_expt();
glmodel = 143; % this is the one that has a regressor for each trial onset
subjs = getGoodSubjects();
runs = 1:metadata.runsPerSubject;

if strcmp(regressor_prefix, 'trial_onset')
    trials = 1:metadata.trialsPerRun;
elseif strcmp(regressor_prefix, 'feedback_onset')
    trials = 1:metadata.trainingTrialsPerRun;
end
    

tmaps = [];

for subj = subjs
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    
    for run = runs
        for trial = trials
            regressor = ['Sn(', num2str(run), ') ', regressor_prefix, '_', num2str(trial)];
            
            tmap = ccnl_get_tmap(EXPT, glmodel, regressor, mask, subj);
            assert(size(tmap, 1) == 1);
            %beta(isnan(beta)) = 0; 
            
            which_row = strcmp(data.participant, metadata.allSubjects{subj}) & data.runId == run & data.newTrialId == trial;
            assert(sum(which_row) == 1);
            
            % initialize betas array if we haven't yet
            %
            if isempty(tmaps)
                tmaps = nan(size(data.participant, 1), size(tmap, 2));
            end
            % append beta
            tmaps(which_row, :) = tmap;
            
            fprintf('loaded subj %d, run %d, trail %d (%s): %d voxels\n', subj, run, trial, regressor, sum(~isnan(tmap)));
        end
    end
end
