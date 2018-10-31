function betas = load_raw_betas(regressor_prefix, data, metadata)

% Load all the betas corresponding to activations at trial onset for the
% given regressor prefix. Does this for all trials for all runs for all subjects.
%
% INPUT:
% regressor_prefix = 'trial_onset' or 'feedback_onset'; assumes regressor
%                    names are of the form {regressor prefix}_{trial id}
% data, metadata = subject behavioral data as output by load_data
%
% OUTPUT:
% betas = [nRows x nVoxels] beta coefficients for each trial.
%         Each row of this matrix corresponds to a single trial == a single
%         row in data.
%

EXPT = context_expt();
glmodel = 143; % this is the one that has a regressor for each trial onset
subjs = getGoodSubjects();
runs = 1:metadata.runsPerSubject;
trials = 1:metadata.trialsPerRun;

betas = [];

for subj = subjs
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));

    for run = runs
        for trial = trials
            regressor = ['Sn(', num2str(run), ') ', regressor_prefix, '_', num2str(trial)];

            beta = [];

            % Find beta for regressor with given name (from ccnl_get_beta)
            %
            for i = 1:length(SPM.xX.name)
                if ~isempty(strfind(SPM.xX.name{i},[regressor,'*'])) || ~isempty(strfind(SPM.xX.name{i},[regressor,'^']))
                    V = spm_vol(fullfile(modeldir,sprintf('beta_%04d.nii',i)));    % residual variance image
                    Y = spm_read_vols(V);
                    
                    beta = Y(:)';
                end
            end            

            assert(~isempty(beta));
            save('test.mat');
            assert(size(beta, 1) == 1);
            %beta(isnan(beta)) = 0;

            which_row = strcmp(data.participant, metadata.allSubjects{subj}) & data.runId == run & data.newTrialId == trial;
            assert(sum(which_row) == 1);

            % initialize betas array if we haven't yet
            %
            if isempty(betas)
                betas = nan(size(data.participant, 1), size(beta, 2));
            end
            % append beta
            betas(which_row, :) = beta;

            fprintf('loaded subj %d, run %d, trail %d (%s): %d voxels\n', subj, run, trial, regressor, sum(~isnan(beta)));
        end
    end
end

