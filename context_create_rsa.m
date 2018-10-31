function rsa = context_create_rsa(rsa_idx, subj)

    % Create rsa structure, helper function for creating EXPT in
    % context_expt.m
    %
    % USAGE: rsa = context_create_rsa(model,subj)
    %
    % INPUTS:
    %   rsa_idx - positive integer indicating which RSA we're doing
    %   subj - integer specifying which subject is being analyzed
    %
    % OUTPUTS:
    %   rsa - a structure with the following fields:
    %     .glmodel - which GLM to use to get the trial-by-trial betas; make sure to have a unique regressor for each trial, e.g. 'trial_onset_1', 'trial_onset_2', etc.
    %     .event - which within-trial event to use for neural activity; used to pick the right betas (needs to be substring of the regressor name), e.g. 'trial_onset'
    %     .mask - path to .nii file, or 3D binary vector of voxels to consider
    %     .radius - searchlight radius in voxels
    %     .which_betas - logical mask for which betas (trials) front the GLM to include (e.g. not timeouts)
    %     .model - struct array describing the models used for behavioral RDMs (see Kriegeskorte et al. 2008) with the fields:
    %         .name - model name
    %         .features - [nTrials x D] feature vector
    %         .distance_measure - name (e.g. 'cosine') or function handler to be used as a distance measure for the RDMs (passed to pdist, see MATLAB documentation)
    %         .is_control - whether this is a control model (e.g. time)
    %
    % Momchil Tomov, Sep 2018


    fprintf('rsa %d, subj %d\n', rsa_idx, subj);

    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
    
    [allSubjects, subjdirs, nRuns] = context_getSubjectsDirsAndRuns();
    assert(isequal(allSubjects, metadata.allSubjects));

    which_trials = data.which_rows & data.isTrain & strcmp(data.participant, allSubjects{subj});; % Look at training trials only
    
    fprintf('which_trials = %s\n', sprintf('%d', which_trials));

    %% Simulate behavior using Kalman filter
    %
    [params, which_structures] = model_params('results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat')
    which_structures = logical(which_structures);
    simulated = simulate_subjects(data, metadata, params, which_structures);

    % RSAs
    %
    switch rsa_idx

        % from the paper model
        %
        case 1
            rsa.event = 'feedback_onset';
            rsa.glmodel = 143;
            rsa.radius = 4 / 1.5;
            rsa.mask = 'masks/mask.nii';
            rsa.which_betas = logical(ones(1, sum(which_trials)));

            rsa.model(1).name = 'posterior';
            rsa.model(1).features = simulated.P(which_trials, which_structures);;
            rsa.model(1).distance_measure = 'cosine';
            rsa.model(1).is_control = false;

            % controls

            trial_onset = cellfun(@str2double, data.actualChoiceOnset);
            rsa.model(2).name = 'time';
            rsa.model(2).features = trial_onset(which_trials);
            rsa.model(2).distance_measure = 'euclidean';
            rsa.model(2).is_control = true;

            rsa.model(3).name = 'run';
            rsa.model(3).features = data.runId(which_trials);
            rsa.model(3).distance_measure = @(c1, c2) c1 ~= c2;
            rsa.model(3).is_control = true;

        otherwise
            assert(false, 'invalid rsa_idx -- should be one of the above');

    end % end of switch statement

end
