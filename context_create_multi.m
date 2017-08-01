function multi = context_create_multi(glmodel, subj, run, save_output)

    % Create multi structure, helper function for creating EXPT in
    % imageryExpt.m
    %
    % USAGE: multi = imagery_create_multi(model,subj,run)
    %
    % INPUTS:
    %   glmodel - positive integer indicating general linear model
    %   subj - integer specifying which subject is being analyzed
    %   run - integer specifying the run
    %
    % OUTPUTS:
    %   multi - a structure with the folloowing fields
    %        .names{i}
    %        .onsets{i}
    %        .durations{i}
    %        optional:
    %        .pmod(i).name
    %        .pmod(i).param
    %        .pmod(i).poly
    %
    % Cody Kommers, July 2016
    
    if nargin < 4 || isempty(save_output)
        save_output = false;
    end
    
    utils; % include some functions
    
    fprintf('glm %d, subj %d, run %d\n', glmodel, subj, run);
    
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
    
    [allSubjects, subjdirs, nRuns] = context_getSubjectsDirsAndRuns();
    assert(isequal(allSubjects, metadata.allSubjects));
        
    % pick the trials that correspond to that subject & run
    % notice that we're not &-ing with data.which_rows -- we might want to
    % run the GLM even for subjects that are not good. In fact, we cannot
    % determine whether subjects are good or not before we have run the GLM
    % and inspected for stuff like rotation and such.
    %
    which_train = ~data.drop & data.isTrain & strcmp(data.participant, allSubjects{subj}) & data.runId == run;
    which_test = ~data.drop & ~data.isTrain & strcmp(data.participant, allSubjects{subj}) & data.runId == run;
    which_rows = which_test | which_train;
    assert(sum(which_train) == metadata.trainingTrialsPerRun);
    assert(sum(which_test) == metadata.testTrialsPerRun);
    assert(sum(which_rows) == metadata.trialsPerRun);
    
    % ...never mind what the thing above said;
    % we only support good subjects here now
    %
    assert(ismember(subj, getGoodSubjects()));
    

    
    % condition = context role for the run
    %
    condition = data.contextRole(which_train);
    condition = condition{1};
    
    % Run model on training trials
    %
    
    if glmodel < 127 || glmodel >= 137
        % Most models use the fixed effects parameter fit with the pilot
        % data
        %
        prior_variance = 0.1249;
        inv_softmax_temp = 2.0064;
        params = [prior_variance inv_softmax_temp];
        which_structures = [1 1 1 0];
    else
        % Only a small subset of the models use the random effects parameter fit on the fMRI data 
        %
        load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts.mat'), 'results', 'results_options');
        params = results(1).x;
        options = results_options(1);
        disp('Using parameters:');
        disp(params);
        disp('generated with options:');
        disp(options);
        assert(options.isFmriData == true);
        assert(~options.fixedEffects);
        assert(isequal(options.which_structures, [1 1 1 0]));
        which_structures = options.which_structures;
    end
    
    simulated = simulate_subjects(data, metadata, params, which_structures, which_rows);        
    
    % get the data for training trials on this run only
    %
    outcomes = data.outcome(which_train);
    choices = simulated.pred(which_train, :);
    P = simulated.P(which_train, :);
    values = simulated.values(which_train, :);
    valuess = simulated.valuess(which_train, :);
    likelihoods = simulated.likelihoods(which_train, :);
    new_values = simulated.new_values(which_train, :);
    new_valuess = simulated.new_valuess(which_train, :);
    lambdas = simulated.lambdas(which_train, :);
    
    % calculate entropy
    % exclude M4 which has P = 0
    %
    H = - sum(P(:,1:3) .* log(P(:, 1:3)), 2);
    H(isnan(H)) = 0; % if a posterior is 0, the entropy is 0
    
    % get the data for test trials on this run only
    %
    test_choices = simulated.pred(which_test, :);
    test_values = simulated.values(which_test, :);
    test_valuess = simulated.valuess(which_test, :);

    % Parametric modulators
    %
    switch glmodel
        
        % cue and outcome regressor for each model
        
        % Simple. Contrast conditions at feedback time
        %
        %
        % neural correlate of latent structure prob
        % show it's predictive of the test performance
        % posterior over structures in OFC and vmPFC
        % hippo encodes context sensitity
        % e.g. hippo only active in modulatory condition -- 
        %     or hippo just codes context, regardless of functional
        %     role
        %     whether there's significant info in hippo re context
        %    option -- classify which context based on hippo activity
        %        separate GLM --
        %           regressor for each context at time of stimulus 
        %           to avoid collinearity, no stimulus event regressor
        %           just events, no pmods
        %           and a feedback event
        %           for each subject, pull out voxel activity in hippo
        %           leave 1 of 9 runs out
        %           gives you 1 beta map for each context per run
        %           have 16 training examples 
        %           train a classifier on 16 examples 
        %           
        %           simple thing -- is context activated selectively in
        %           modulatory condition
        %               vs. irrelevant and additive
        %           modulatory - irrelevant
        %           additive - irrelevant
        %           modulatory - additive
        %           
        %
        % compare 1 (feedback) with 140 (trial onset) with 141 (null)
        %
        case 1 % <------------- GOOD
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M2 posterior pmod @ outcome + extra stuff for test trials ZOMG
        %
        case 2 % <------------- GOOD
            % M2 (modulatory) posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_posterior';
            multi.pmod(1).param{1} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % Events 3, 4, 5, 6: Test stimulus, one per
            % condition-cue-context tuple
            %
            % name = test_{condition}_x{data.cueId + 1}c{data.contextId + 1}
            %   e.g. test_irrelevant_x1c3
            % onset = stimulus onset
            % duration = 0
            %
            % identify voxels that predict expected outcome purely on training
            % on test trials, at stimulus onset, do they show the pattern
            % of modulation that we see in behavior
            % just look at the betas then (no pmods for test trials)
            %
            test_data.cueIds = data.cueId(which_test);
            test_data.contextIds = data.contextId(which_test);
            test_data.actualChoiceOnsets = cellfun(@str2num, data.actualChoiceOnset(which_test))';
            for i=1:length(test_data.actualChoiceOnsets)
                multi.names{2 + i} = sprintf('test_%s_x%dc%d', condition, test_data.cueIds(i) + 1, test_data.contextIds(i) + 1);
                multi.onsets{2 + i} = test_data.actualChoiceOnsets(i);
                multi.durations{2 + i} = 0;
            end

        % correct vs. wrong pmod @ outcome
        % TODO WTFFFFFFFFFFFFFFFFFFFFFFFFFF didn't work on NCF
        %
        case 3
            % correct vs. wrong (1/0) @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'correct';
            multi.pmod(1).param{1} = data.response.corr(which_train)'; % correct (1 or 0) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % prediction error pmod @ outcome
        % WRONG! choices != expected outcome
        %
        case 4
            % prediction error @ feedback onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - choices'; % outcome - expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % outcome (sick vs. not sick) @ outcome
        %
        case 5 % GOOD
            % sick vs. not sick @ feedback / outcome onset (trials 1..20)
            %
            multi.names{1} = 'sick';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train & strcmp(data.sick, 'Yes')))';
            multi.durations{1} = zeros(size(data.contextRole(which_train & strcmp(data.sick, 'Yes'))));            
            
            multi.names{2} = 'not sick';
            multi.onsets{2} = cellfun(@str2num, data.actualFeedbackOnset(which_train & ~strcmp(data.sick, 'Yes')))';
            multi.durations{2} = zeros(size(data.contextRole(which_train & ~strcmp(data.sick, 'Yes'))));                        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));
            
        % outcome pmod @ outcome -- sanity check
        % should be same as 5
        % result: nothing after FWE; asymmetrical V1
        %
        case 6 % GOOD
            % sick vs. not sick @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            assert(isequal(strcmp(data.sick(which_train), 'Yes'), outcomes));
            multi.pmod(1).name{1} = 'outcome';
            multi.pmod(1).param{1} = outcomes'; % outcome == sick (1 or 0) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % expected outcome pmod @ start
        % WRONG! choices != expected outcome
        %
        case 7
            % expected outcome @ trial onset (trials 1..20)
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'expected_outcome';
            multi.pmod(1).param{1} = choices'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % const @ feedback onset
            %
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M3 posterior pmod @ outcome
        %
        case 8 % GOOD
            % M3 (additive) posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_posterior';
            multi.pmod(1).param{1} = P(:,3)'; % posterior P(M3 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M1 posterior pmod @ outcome
        %
        case 9 % GOOD
            % M1 (irrelevant) posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_posterior';
            multi.pmod(1).param{1} = P(:,1)'; % posterior P(M1 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % posterior entropy @ outcome
        %
        case 10 % OKAY
            % posterior entropy @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'posterior_entropy';
            multi.pmod(1).param{1} = H'; % entropy of P(M | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % posterior entropy + context role / condition @ outcome
        %
        case 11
            % posterior entropy @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = strcat(condition, '_posterior_entropy');
            multi.pmod(1).param{1} = H'; % entropy of P(M | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % flipped posterior entropy @ outcome
        %
        case 12
            % flipped posterior entropy @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'flipped_posterior_entropy';
            multi.pmod(1).param{1} = max(H') - H'; % flipped entropy of P(M | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % flipped posterior entropy + context role / condition @ outcome
        %
        case 13
            % flipped posterior entropy @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = strcat(condition, '_flipped_posterior_entropy');
            multi.pmod(1).param{1} = max(H') - H'; % flipped entropy of P(M | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % posterior std @ outcome
        %
        case 14
            % posterior std @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'posterior_std';
            multi.pmod(1).param{1} = std(P(:, 1:3)'); % std of P(M | h_1:n) for trials 1..20, exluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % posterior std + context role / condition @ outcome
        %
        case 15
            % posterior std @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = strcat(condition, '_posterior_std');
            multi.pmod(1).param{1} = std(P(:, 1:3)'); % std of P(M | h_1:n) for trials 1..20, exluding M4
            multi.pmod(1).poly{1} = 1; % first order        
            
            % expected outcome @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M2 posterior - M1 posterior pmod @ outcome
        %
        case 16
            % M2 - M1 posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_minus_M1_posterior';
            multi.pmod(1).param{1} = P(:,2)' - P(:,1)';
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M3 posterior - M1 posterior pmod @ outcome
        %
        case 17
            % M3 - M1 posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_minus_M1_posterior';
            multi.pmod(1).param{1} = P(:,3)' - P(:,1)';
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M2 & M1 posterior pmods @ outcome
        % result: nothing BUT also do the individual regressors! omg...
        %
        case 18
            % M2 (modulatory) & M1 (irrelevant) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M2_posterior';
            multi.pmod(1).param{1} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'M1_posterior';
            multi.pmod(1).param{2} = P(:,1)'; % posterior of P(M1 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M3 & M1 posterior pmods @ outcome
        %
        case 19
            % M3 (additive) & M1 (irrelevant) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M3_posterior';
            multi.pmod(1).param{1} = P(:,3)'; % posterior P(M3 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'M1_posterior';
            multi.pmod(1).param{2} = P(:,1)'; % entropy of P(M1 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M3 & M2 posterior pmods @ outcome
        %
        case 20
            % M3 (additive) & M2 (modulatory) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M3_posterior';
            multi.pmod(1).param{1} = P(:,3)'; % posterior P(M3 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'M2_posterior';
            multi.pmod(1).param{2} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % Response @ reaction onset 
        %
        case 21
            % button press @ reaction onset (trials 1..20)
            % 
            multi.names{1} = 'keypress';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOffset(which_train & ~data.no_response))';
            multi.durations{1} = zeros(size(data.contextRole(which_train & ~data.no_response)));
            
            multi.pmod(1).name{1} = 'pressed_sick';
            multi.pmod(1).param{1} = strcmp(data.response.keys(which_train & ~data.no_response), 'left')'; % whether subject pressed "sick", for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ outcome / feedback onset (trials 1..20)
            % 
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train & ~data.no_response))';
            multi.durations{2} = zeros(size(data.contextRole(which_train & ~data.no_response)));
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train & ~data.no_response))';
            multi.durations{3} = zeros(size(data.contextRole(which_train & ~data.no_response)));
       
            
        % Blocked M2 & M1 posterior pmods @ outcome
        %
        case 22
            % M2 (modulatory) & M1 (irrelevant) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = strcat('M2_posterior_', condition);
            multi.pmod(1).param{1} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = strcat('M1_posterior_', condition);
            multi.pmod(1).param{2} = P(:,1)'; % entropy of P(M1 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % Blocked M3 & M1 posterior pmods @ outcome
        %
        case 23
            % M3 (additive) & M1 (irrelevant) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = strcat('M3_posterior_', condition);
            multi.pmod(1).param{1} = P(:,3)'; % posterior P(M3 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = strcat('M1_posterior_', condition);
            multi.pmod(1).param{2} = P(:,1)'; % entropy of P(M1 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % Blocked M3 & M2 posterior pmods @ outcome
        %
        case 24
            % M3 (additive) & M2 (modulatory) posteriors @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = strcat('M3_posterior_', condition);
            multi.pmod(1).param{1} = P(:,3)'; % posterior P(M3 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = strcat('M2_posterior_', condition);
            multi.pmod(1).param{2} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20, excluding M4
            multi.pmod(1).poly{2} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % Expected outcome @ RT
        % Actual outcome @ feedback
        % result: expected -- strange bilateral in ventricles (behind
        %                     hippo), after cluster FWE
        %         actual -- lots: negative of above, bilateral insula,
        %                   bilateral OFC?
        %
        case 25
            % Expected outcome @ RT
            %
            multi.names{1} = 'expected_outcome';
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num,RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'expected';
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % Actual outcome @ feedback
            %
            multi.names{2} = 'actual_outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'actual';
            multi.pmod(2).param{1} = outcomes'; % outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order  
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));
            
        % Prediction error @ feedback
        % result: all over the place but not where we want it
        %
        case 26
            % Prediction error @ feedback
            %
            multi.names{1} = 'prediction_error';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'actual';
            multi.pmod(1).param{1} = outcomes' - values'; % outcome - expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M2 posterior pmod @ outcome
        % result: a bunch of activation all over
        %         the problem is, 'M2_posterior' very strongly correlates with
        %         'feedback' and is mostly orthogonalized out
        %
        case 27 % <------------- GOOD
            % M2 (modulatory) posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_posterior';
            multi.pmod(1).param{1} = P(:,2)'; % posterior P(M2 | h_1:n) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            

        % context role @ trial onset
        % WRONG! cannot have 2 params with the same onset!
        %
        case 28
            % context role @ trial onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % context role @ trial onset, no separate trial_onset regressor
        % (almost same as 28)
        % result: don't get hippocampus in 'additive - irrelevant' !
        %
        case 29
            % context role @ trial onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));            
            
            
        % --------------- PREDICTED VALUES, multiple models ------------
            
            
        % M1, M2 & M3 value pmod s @ trial onset (before updated)
        % WRONG -- duplicate onsets...
        %
        case 30
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_value';
            multi.pmod(1).param{2} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_value';
            multi.pmod(1).param{3} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M1, M2 & M3 value pmods @ trial onset
        % without the separate trial_onset regressor
        % (almost same as 30)
        %
        case 31
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_value';
            multi.pmod(1).param{2} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_value';
            multi.pmod(1).param{3} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        

        % M1, M2 & M3 value pmods @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % WRONG -- duplicate onsets...
        %
        case 32
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_value';
            multi.pmod(1).param{2} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_value';
            multi.pmod(1).param{3} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - values'; % PE = outcome - expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));

        % M1, M2 & M3 value pmods @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 32; also look a 30/31)
        % result: prediction_error -- all over the place; big blob in mPFC
        %
        case 33
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_value';
            multi.pmod(1).param{2} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_value';
            multi.pmod(1).param{3} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - values'; % PE = outcome - expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      
            
         
        % ---------------- PREDICTED VALUES + PREDICTION ERRORS, single model ----------------
            
            
        % M1 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % WRONG -- duplicate onsets...
        %
        case 34
            % M1 (irrelevant) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,1)'; % PE = outcome - expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));

        % M2 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % WRONG -- duplicate onsets...
        %
        case 35
            % M2 (modulatory) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_value';
            multi.pmod(1).param{1} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,2)'; % PE = outcome - expected outcome by M2 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));

        % M3 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % WRONG -- duplicate onsets...
        %
        case 36
            % M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_value';
            multi.pmod(1).param{1} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,3)'; % PE = outcome - expected outcome by M3 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));

            
        % M1 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 34)
        % result: nothing really
        %
        case 37
            % M1 (irrelevant) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,1)'; % PE = outcome - expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      

        % M2 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 35)
        % result: b/4 FWE -- 'M2_value' -> weird visual neg
        %
        case 38
            % M2 (modulatory) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_value';
            multi.pmod(1).param{1} = valuess(:,2)'; % values predicted by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,2)'; % PE = outcome - expected outcome by M2 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order      

        % M3 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 36)
        % result: b/4 FWE -- 'M3_value' -> bilateral posterior hippocampus
        %                                  weird visual neg
        %
        case 39
            % M3 (additive) predicted values @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_value';
            multi.pmod(1).param{1} = valuess(:,3)'; % values predicted by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % Prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,3)'; % PE = outcome - expected outcome by M3 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order   
            
            
        % --------------- LIKELIHOODS ---------------------
            
        
        % M1, M2 & M3 likelihood pmods @ feedback (outcome) onset (before updated)
        %
        case 40
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) outcome likelihood @ feedback (outcome) onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_lik';
            multi.pmod(1).param{1} = likelihoods(:,1)'; % outcome likelihoods predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_lik';
            multi.pmod(1).param{2} = likelihoods(:,2)'; % outcome likelihoods predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_lik';
            multi.pmod(1).param{3} = likelihoods(:,3)'; % outcome likelihoods predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M1, M2 & M3 likelihood pmods @ feedback (outcome) onset (before updated)
        % without the separate trial_onset regressor
        % (almost same as 40)
        %
        case 41
            % M1 (irrelevant), M2 (modulatory) & M3 (additive) outcome likelihood @ feedback (outcome) onset (trials 1..20)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.orth{1} = 0; % do NOT orthogonalize them!
            
            multi.pmod(1).name{1} = 'M1_lik';
            multi.pmod(1).param{1} = likelihoods(:,1)'; % outcome likelihoods predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            multi.pmod(1).name{2} = 'M2_lik';
            multi.pmod(1).param{2} = likelihoods(:,2)'; % outcome likelihoods predicted by M2 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            multi.pmod(1).name{3} = 'M3_lik';
            multi.pmod(1).param{3} = likelihoods(:,3)'; % outcome likelihoods predicted by M3 for trials 1..20
            multi.pmod(1).poly{3} = 1; % first order        
            
        
        % --------------- PREDICTION ERRORS, multi-model ----------------
        

        % PE @ feedback (outcome) onset
        % (also look at 32, 30/31)
        %
        case 42
            % Prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - values'; % PE = outcome - expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 42; also look a 33, 32, 30/31)
        %
        case 43
            % Prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - values'; % PE = outcome - expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      
            
            
                
        % ----------- PREDICTION ERRORS, single model ---------------    
           
        
        
        % M1 only PE @ feedback (outcome) onset
        %
        case 44
            % M1 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,1)'; % PE = outcome - expected outcome by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M2 only PE @ feedback (outcome) onset
        %
        case 45
            % M2 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,2)'; % PE = outcome - expected outcome by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M3 only PE @ feedback (outcome) onset
        %
        case 46
            % M3 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,3)'; % PE = outcome - expected outcome by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % M1 only PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 44)
        %
        case 47
            % M1 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,1)'; % PE = outcome - expected outcome by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      

        % M2 only PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 45)
        %
        case 48
            % M2 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,2)'; % PE = outcome - expected outcome by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      

        % M3 only PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 46)
        %
        case 49
            % M3 prediction error @ feedback
            %
            multi.names{1} = 'outcome';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'prediction_error';
            multi.pmod(1).param{1} = outcomes' - valuess(:,3)'; % PE = outcome - expected outcome by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order      
            
        % Value (expected outcome) @ trial onset
        % Prediction error @ feedback (outcome) onset
        % WRONG -- 'trial_onset' duplicates 'stimulus'
        %
        case 50
            % Value (expected outcome) @ trial onset
            %
            multi.names{1} = 'stimulus';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'value';
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % Prediction error @ feedback (outcome) onset
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - values'; % outcome - expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order  
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));
            
        % Value (expected outcome) @ trial onset
        % Prediction error @ feedback (outcome) onset
        % w/o the visual stimulus regressor
        % (almost same as 50)
        % result: 'value' => nothing, 'prediction_error' => all over
        %
        case 51
            % Value (expected outcome) @ trial onset
            %
            multi.names{1} = 'stimulus';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'value';
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % Prediction error @ feedback (outcome) onset
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - values'; % outcome - expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order  
            
        % M3 posterior - M2 posterior pmod @ outcome
        %
        case 52
            % M3 - M2 posterior @ feedback / outcome onset (trials 1..20)
            % 
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_minus_M2_posterior';
            multi.pmod(1).param{1} = P(:,3)' - P(:,2)';
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % Bayesian surprise = Kullback?Leibler divergence @ feedback
        % (outcome) onset
        % https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
        % result: bilteral rlPFC, after cluster FWE
        %         negative -- mPFC, posterior cingulate
        %
        case 53
            % Q(t,:) = prior distribution over causal structures (models)
            % at trial t (i.e. before update)
            % P(t,:) = posterior (i.e. after update)
            % Note that Q(t+1,:) = P(t,:)
            % surprise(t) = D_KL(P(t,:) || Q(t,:)) = sum P(t,i) log( P(t,i) / Q(t,i) )
            % 
            priors = which_structures / sum(which_structures);
            Q = [priors; P(1:end-1,:)];
            logs = log2(P) - log2(Q); 
            logs(isnan(logs)) = 0; % lim_{x->0} x log(x) = 0
            surprise = sum(P .* logs, 2);
            surprise(isnan(surprise)) = 0; % weird things happen when P --> 0, e.g. we get -Infs

            % sanity check
            %logs2 = log2(P ./ Q); 
            %logs2(isnan(logs2)) = 0; % lim_{x->0} x log(x) = 0
            %disp(P)
            %disp(logs)
            %disp(logs2)
            %logs2 - logs
            %assert(abs(sum(sum(logs2(~isinf(logs2)) - logs(~isinf(logs))))) < 1e-6);
            
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'surprise';
            multi.pmod(1).param{1} = surprise';
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % One regressor for each training trial onset (for classifier)
        % http://ufldl.stanford.edu/tutorial/supervised/ExerciseConvolutionalNeuralNetwork/
        %
        case 54
            trial_onsets = data.actualChoiceOnset(which_train);
            for t=1:20
                multi.names{t} = ['trial_onset_', num2str(t)];
                multi.onsets{t} = [str2double(trial_onsets(t))];
                multi.durations{t} = [0];
            end
            
        % One regressor for each test trial onset (for classifier)
        % http://ufldl.stanford.edu/tutorial/supervised/ExerciseConvolutionalNeuralNetwork/
        %
        case 55
            trial_onsets = data.actualChoiceOnset(which_test);
            for t=1:4
                multi.names{t} = ['trial_onset_', num2str(t)];
                multi.onsets{t} = [str2double(trial_onsets(t))];
                multi.durations{t} = [0];
            end            
            
        % main effect but 1 s durations
        % result: additive - irrelevant becomes less significant
        case 56
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train)));

        % main effect but 3 s durations at trial onset
        % result: nothing for additive - irrelevant
        case 57
            % context role @ trial onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = 3 * ones(size(data.contextRole(which_train)));            

        % main effect but 1 s durations at ITI (sanity check)
        % result: nothing for additive - irrelevant
        case 58
            % context role @ trial onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualItiOnset(which_train))';
            multi.durations{1} = 1 * ones(size(data.contextRole(which_train)));            
            
        % One regressor for each training feedback onset (for classifier)
        % http://ufldl.stanford.edu/tutorial/supervised/ExerciseConvolutionalNeuralNetwork/
        %
        case 59
            feedback_onsets = data.actualFeedbackOnset(which_train);
            for t=1:20
                multi.names{t} = ['trial_onset_', num2str(t)];
                multi.onsets{t} = [str2double(feedback_onsets(t))];
                multi.durations{t} = [0];
            end

        % One regressor for each training + test trial onset (for classifier)
        % same as 54 + 55
        % http://ufldl.stanford.edu/tutorial/supervised/ExerciseConvolutionalNeuralNetwork/
        %
        case 60
            % training trials
            trial_onsets = data.actualChoiceOnset(which_train);
            for t=1:20
                multi.names{t} = ['trial_onset_', num2str(t)];
                multi.onsets{t} = [str2double(trial_onsets(t))];
                multi.durations{t} = [0];
            end
            
            % test trials
            trial_onsets = data.actualChoiceOnset(which_test);
            for t=21:24
                multi.names{t} = ['trial_onset_', num2str(t)];
                multi.onsets{t} = [str2double(trial_onsets(t - 20))];
                multi.durations{t} = [0];
            end            

            
        % Value (expected outcome) @ trial onset, separated by condition
        % result: before FWE, 'value_additive - value_irrelevant' -> bilateral putamen
        %                     'value_additive - value_modulatory' -> putamen, hippocampus (flip signs)
        %
        case 61
            % Value (expected outcome) @ trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Value (expected outcome) @ RT, separated by condition
        % result: before FWE, 'value_modulatory - value_irrelevant' -> S1
        %                     'value_additive - value_irrelevant' -> R Put
        %
        case 62
            % Value (expected outcome) @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % Cause-specific value @ trial onset, separated by condition
        % result: nothing
        %
        case 63
            % Cause-specific value (expected outcome) @ trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            multi.pmod(1).name{1} = ['M', num2str(M), '_value_', condition];
            multi.pmod(1).param{1} = valuess(:,M)'; % cause-specific value (expected outcome) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Cause-specific value @ RT, separated by condition
        % result: nothing
        %
        case 64
            % Cause-specific value (expected outcome) @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            multi.pmod(1).name{1} = ['M', num2str(M), '_value_', condition];
            multi.pmod(1).param{1} = valuess(:,M)'; % cause-specific value (expected outcome) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % main effect @ feedback for trials 10..20
        % result: almost same as main effect on 1..20
        %
        case 65
            which_trials = which_train & data.trialId >= 10;
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));
            
        % main effect @ RT, trials 10..20
        % result: almost same as main effect on 1..20
        %
        case 66
            which_trials = which_train & data.trialId >= 10;
            
            % context role @ RT
            % 
            multi.names{1} = condition;
            RTs = data.actualChoiceOffset(which_trials);
            feedback_time = data.actualFeedbackOnset(which_trials);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));

        % Response @ 1 s before reaction onset 
        % result: same as pressed_sick @ RT => big opposite blobs in visual
        %
        case 67
            % button press @ 1 s before reaction onset (trials 1..20)
            % 
            multi.names{1} = 'keypress';
            multi.onsets{1} = -1 + cellfun(@str2num,data.actualChoiceOffset(which_train & ~data.no_response))';
            multi.durations{1} = zeros(size(data.contextRole(which_train & ~data.no_response)));
            
            multi.pmod(1).name{1} = 'pressed_sick';
            multi.pmod(1).param{1} = strcmp(data.response.keys(which_train & ~data.no_response), 'left')'; % whether subject pressed "sick", for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ outcome / feedback onset (trials 1..20)
            % 
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train & ~data.no_response))';
            multi.durations{2} = zeros(size(data.contextRole(which_train & ~data.no_response)));
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train & ~data.no_response))';
            multi.durations{3} = zeros(size(data.contextRole(which_train & ~data.no_response)));

        % Response @ trial onset 
        % result: same as pressed_sick @ RT => big opposite blobs in visual
        %
        case 68
            % button press @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train & ~data.no_response))';
            multi.durations{1} = zeros(size(data.contextRole(which_train & ~data.no_response)));
            
            multi.pmod(1).name{1} = 'pressed_sick';
            multi.pmod(1).param{1} = strcmp(data.response.keys(which_train & ~data.no_response), 'left')'; % whether subject pressed "sick", for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ outcome / feedback onset (trials 1..20)
            % 
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train & ~data.no_response))';
            multi.durations{2} = zeros(size(data.contextRole(which_train & ~data.no_response)));

        % Value (expected outcome) @ 1 s after trial onset, separated by condition
        % result: before FWE, 'value_modulatory - value_irrelevant' -> S1 
        %         before FWE, 'value_modulatory - value_additive' -> left
        %                      hippocampus, P = 0.01
        %
        case 69
            % Value (expected outcome) @ 1 s after trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Value (expected outcome) @ 1 s before RT, separated by condition
        % result: before FWE, 'value_additive - value_irrelevant' --> bilateral putamen, some hippo
        %
        case 70
            % Value (expected outcome) 1 s before @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        %
        % ----------- THESE INCLUDE THE TEST TRIALS -----------------
        %
            
        % main effect for ALL trials (including test)
        % result: before FWE, 'modulatory - irrelevant' -> parietal
        %                     'modulatory - additive' -> putamen (neg)
        %                     'additive - irrelevant' -> midbrain?     
        %
        case 71
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
                    
        % Value (expected outcome) @ trial onset, separated by condition INCLUDING test trials
        % same as 61
        % result: similar to 61: 
        %         before FWE, 'value_modulatory - value_additive' -> putamen, hippocampus (flip signs)
        %
        case 72
            % Value (expected outcome) @ trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = [values' test_values']; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Value (expected outcome) @ RT, separated by condition INCLUDING test trials
        % same as 62
        % result: nothing
        %
        case 73
            % Value (expected outcome) @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = [values' test_values']; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
        % Cause-specific value @ trial onset, separated by condition INCLUDING test trials
        % same as 63
        % result: stupid ; this is the same as just value
        %
        case 74
            % Cause-specific value (expected outcome) @ trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            multi.pmod(1).name{1} = ['M', num2str(M), '_value_', condition];
            multi.pmod(1).param{1} = [valuess(:,M)' test_valuess(:,M)']; % cause-specific value (expected outcome) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Cause-specific value @ RT, separated by condition INCLUDING test trials
        % same as 64
        % result: stupid ; this is the same as just value
        %
        case 75
            % Cause-specific value (expected outcome) @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            multi.pmod(1).name{1} = ['M', num2str(M), '_value_', condition];
            multi.pmod(1).param{1} = [valuess(:,M)' test_valuess(:,M)']; % cause-specific value (expected outcome) for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
        % main effect @ feedback for all trials INCLUDING test trials
        % result: same as 1. Nothing really
        %
        case 76
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
        % main effect @ RT, all trials INCLUDING test trials
        % result: nothing interesting
        %
        case 77
            % context role @ RT
            % 
            multi.names{1} = condition;
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));

        % Response @ 1 s before reaction onset  INCLUDING test trials
        % same as 67
        % result: big visual activations corresponding to where the sick
        %         face is on the screen
        %
        case 78
            % button press @ 1 s before reaction onset (trials 1..20)
            % 
            multi.names{1} = 'keypress';
            multi.onsets{1} = -1 + cellfun(@str2num,data.actualChoiceOffset(which_rows & ~data.no_response))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows & ~data.no_response)));
            
            multi.pmod(1).name{1} = 'pressed_sick';
            multi.pmod(1).param{1} = strcmp(data.response.keys(which_rows & ~data.no_response), 'left')'; % whether subject pressed "sick", for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ outcome / feedback onset (trials 1..20)
            % 
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_rows & ~data.no_response))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows & ~data.no_response)));
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_rows & ~data.no_response))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows & ~data.no_response)));

        % Response @ trial onset  INCLUDING test trials
        % same as 68
        % result: similar visual activations as 79 but the positive one
        %         disappeared
        %
        case 79
            % button press @ trial onset (trials 1..20)
            % 
            multi.names{1} = 'keypress';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_rows & ~data.no_response))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows & ~data.no_response)));
            
            multi.pmod(1).name{1} = 'pressed_sick';
            multi.pmod(1).param{1} = strcmp(data.response.keys(which_rows & ~data.no_response), 'left')'; % whether subject pressed "sick", for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ outcome / feedback onset (trials 1..20)
            % 
            multi.names{2} = 'feedback';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_rows & ~data.no_response))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows & ~data.no_response)));
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_rows & ~data.no_response))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows & ~data.no_response)));

        % Value (expected outcome) @ 1 s after trial onset, separated by condition INCLUDING test trials
        % same as 69
        % result: nothing
        %
        case 80
            % Value (expected outcome) @ 1 s after trial onset
            %
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = [values' test_values']; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

        % Value (expected outcome) @ 1 s before RT, separated by condition INCLUDING test trials
        % same as 70
        % result: b/4 FWE, 'value_additive - value_irrelevant' -> bilateral putamen
        %
        case 81
            % Value (expected outcome) 1 s before @ RT
            %
            multi.names{1} = 'RT';
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            multi.pmod(1).name{1} = ['value_', condition];
            multi.pmod(1).param{1} = [values' test_values']; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            

        
        %
        % ------------- MORE MAIN EFFECTS ---------------
        %
        
        
        
        % main effect @ trial onset + 0.5 s, trials 1..20
        % result: nothing
        %
        case 82
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
        % main effect @ trial onset + 1 s, trials 1..20
        % result: nothing
        %
        case 83
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
        % main effect @ RT - 0.5s, trials 1..20
        % result: nothing 
        %
        case 84
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);

            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -0.5 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

        % main effect @ RT - 1s, trials 1..20
        % result: nothing
        %
        case 85
            RTs = data.actualChoiceOffset(which_train);
            feedback_time = data.actualFeedbackOnset(which_train);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);

            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            
            
            
        % main effect @ trial onset + 0.5 s, trials 10..20
        % result: nothing
        %
        case 86
            which_trials = which_train & data.trialId >= 10;
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
        % main effect @ trial onset + 1 s, trials 10..20
        % result: nothing
        %
        case 87
            which_trials = which_train & data.trialId >= 10;
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
        % main effect @ RT - 0.5s, trials 10..20
        % result: nothing
        %
        case 88
            which_trials = which_train & data.trialId >= 10;

            RTs = data.actualChoiceOffset(which_trials);
            feedback_time = data.actualFeedbackOnset(which_trials);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -0.5 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));

        % main effect @ RT - 1s, trials 10..20
        % result: nothing
        %
        case 89
            which_trials = which_train & data.trialId >= 10;
            
            RTs = data.actualChoiceOffset(which_trials);
            feedback_time = data.actualFeedbackOnset(which_trials);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));
            
            
        % main effect @ trial onset + 0.5 s, trials 1..24
        % result: nothing
        %
        case 90
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
        % main effect @ trial onset + 1 s, trials 1..24
        % result: nothing
        %
        case 91
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
        % main effect @ RT - 0.5s, trials 1..24
        % result: nothing
        %
        case 92
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -0.5 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));

        % main effect @ RT - 1s, trials 1..24
        % result: nothing
        %
        case 93
            RTs = data.actualChoiceOffset(which_rows);
            feedback_time = data.actualFeedbackOnset(which_rows);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
        


        % main effect @ trial onset + 0.5 s, trials 10..24
        %
        case 94
            which_trials = (which_train & data.trialId >= 10) | which_test;
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
        % main effect @ trial onset + 1 s, trials 10..24
        %
        case 95
            which_trials = (which_train & data.trialId >= 10) | which_test;
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
        % main effect @ feedback - 1.5s (RT - 0.5s), trials 10..24
        %
        case 96
            which_trials = (which_train & data.trialId >= 10) | which_test;

            RTs = data.actualChoiceOffset(which_trials);
            feedback_time = data.actualFeedbackOnset(which_trials);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -0.5 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));

        % main effect @ feedback - 2s (RT - 1s), trials 10..24
        %
        case 97
            which_trials = (which_train & data.trialId >= 10) | which_test;
            
            RTs = data.actualChoiceOffset(which_trials);
            feedback_time = data.actualFeedbackOnset(which_trials);
            % for the TIMEOUTs, use the time 1 second before feedback (i.e.
            % ISI onset)
            %
            RTs(strcmp(RTs, '')) = cellfun(@num2str, num2cell(cellfun(@str2num, feedback_time(strcmp(RTs, ''))) - 1), 'uniformoutput', 0);

            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = -1 + cellfun(@str2num, RTs)';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));
            
            
        %
        % ------------------ VALUES --------------------
        %
            
        % M4 value pmod @ trial onset (before updated)
        % + PE @ feedback (outcome) onset
        % without the separate trial_onset regressor
        % (almost same as 34)
        % result: before FWE, 'M4_value' -> bilateral putamen
        %         before FWE, 'prediction_error' -> insula, PMA, frontal, negative posterior hippocampus?
        %
        case 98
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M4_value';
            multi.pmod(1).param{1} = valuess(:,4)'; % values predicted by M4 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % M4 prediction error @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(2).name{1} = 'prediction_error';
            multi.pmod(2).param{1} = outcomes' - valuess(:,4)'; % PE = outcome - expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order

        % M1 value + M4 value pmod @ trial onset (before updated)
        % without the separate trial_onset regressor
        % result: cluster FWE, 'M4_value' -> posterior hippo? or ventricles... ??
        %         'M1_value' -> nothing
        %
        case 99
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M4_value';
            multi.pmod(1).param{1} = valuess(:,4)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
                        
            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % values predicted by M4 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order

        % M4 value + M1 value pmod @ trial onset (before updated)
        % without the separate trial_onset regressor
        % same as 99 but flipped, for sanity check
        % result: same as 99
        %
        case 100
            multi.names{1} = 'trial';
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_value';
            multi.pmod(1).param{1} = valuess(:,1)'; % values predicted by M1 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
                        
            multi.pmod(1).name{2} = 'M4_value';
            multi.pmod(1).param{2} = valuess(:,4)'; % values predicted by M4 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
        % learned association with c1, blocked by reported value of c1
        % (x3c1)
        % hypothesis: if you had stronger activation in hippocampus during encoding, your contextual association would be stronger
        %             btw this applies to irrelevant and modulatory only; additive is x1c1+, x2c1+ 
        % expect: 'sick x c1_outcome': significant in hippocampus
        %            i.e. if you thought c1 makes you sick, then you
        %            probably had a stronger activation during the c1 ->
        %            sick trials than the c1 -> not sick trials
        %         'not sick x c1_outcome': not significant in hippocampus
        %            i.e. if you thought c1 doesn't make you sick, you
        %            probably didn't have that => this relies on the
        %            assymetry that people have a prior that things don't
        %            make you sick
        %            OR alternatively, even better, we get the NEGATIVE
        %            activation as sick x c1_outcome!
        %                 ... but then it's confounded with perceiving the actual
        %                 outcome...but I guess how do you disengantle
        %                 encoding & perception?
        %
        %         UGH flipped the sign of 'not_sick x c1_outcome' ...FUCK
        %         so 'c1_outcome' contrast = sick x (encoded sick > not
        %         sick) - not_sick x (encoded not sick > sick)
        %              or alterantive is true => 
        %         but their difference (-) is the actual sum... so if
        %         alternative is true, difference would be very significant
        %
        % result: b/4 FWE, 'sickxc1_outcome - not_sickxc1_outcome' --> left, hippocampus, P = 0.01
        %c
        case 101
            % whether subject thought c1 makes you sick (unlikely but does
            % happen ~25% of the time)
            
            x3c1_choice = data.response.keys(which_test & data.cueId == 2 & data.contextId == 0);
            if strcmp(x3c1_choice, 'left')
                x3c1_choice = 'sick';
            else
                x3c1_choice = 'not_sick';
            end

            which_trials = which_train & data.contextId == 0;
            
            % c1 association @ feedback/outcome onset
            % 
            multi.names{1} = x3c1_choice;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_trials))';
            multi.durations{1} = ones(size(data.contextRole(which_trials))); % 1s durations
            
            % whether c1 made you sick on a given trial
            c1_outcomes = strcmp(data.corrAns(which_trials & data.contextId == 0), 'left');
            if ~strcmp(condition, 'additive') % we don't like additive here (c1 always makes you sick)
                multi.pmod(1).name{1} = 'c1_outcome';
                multi.pmod(1).param{1} = c1_outcomes'; % whether c1 made you sick on trials 1..20
                multi.pmod(1).poly{1} = 1; % first order            
            else
                multi.pmod(1).name{1} = 'c1_outcome_additive';
                multi.pmod(1).param{1} = c1_outcomes'; % whether c1 made you sick on trials 1..20
                multi.pmod(1).poly{1} = 1; % first order            
            end 
            
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));

            
        % learned association with x1, blocked by reported value of x1
        % (x1c3)
        % same as 101 except for x1 instead of c1
        %
        % result: not same activation as 101 => good! hippo is context-dep
        % bla sanity check
        %
        case 102
            % whether subject thought x1 makes you sick (unlikely but does
            % happen ~25% of the time)
            x1c3_choice = data.response.keys(which_test & data.cueId == 0 & data.contextId == 2);
            if strcmp(x1c3_choice, 'left')
                x1c3_choice = 'sick';
            else
                x1c3_choice = 'not_sick';
            end

            which_trials = which_train & data.cueId == 0;
            
            % x1 association @ feedback/outcome onset
            % 
            multi.names{1} = x1c3_choice;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_trials))';
            multi.durations{1} = ones(size(data.contextRole(which_trials))); % 1s durations
            
            % whether c1 made you sick on a given trial
            x1_outcomes = strcmp(data.corrAns(which_trials & data.cueId == 0), 'left');
            if ~strcmp(condition, 'irrelevant') % we don't like irrelevant here (x1 always makes you sick)
                multi.pmod(1).name{1} = 'x1_outcome';
                multi.pmod(1).param{1} = x1_outcomes'; % whether x1 made you sick on trials 1..20
                multi.pmod(1).poly{1} = 1; % first order            
            else
                multi.pmod(1).name{1} = 'x1_outcome_irrelevant';
                multi.pmod(1).param{1} = x1_outcomes'; % whether x1 made you sick on trials 1..20
                multi.pmod(1).poly{1} = 1; % first order            
            end 
           
            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));

       
        % Bayesian surprise = Kullback?Leibler divergence @ feedback
        % (outcome) onset
        % same as 53 but only trials 6..20
        % WRONG: can't do it; not enough info ... contrasts fail sp_conman,
        %                     sp_opp : !invalid contrast
        %
        case 103
            which_trials = which_train & data.trialId >= 6;
            
            priors = which_structures / sum(which_structures);
            Q = [priors; P(1:end-1,:)];
            logs = log2(P ./ Q); 
            logs(isnan(logs)) = 0; % lim_{x->0} x log(x) = 0
            surprise = sum(P .* logs, 2);
            surprise(isnan(surprise)) = 0; % weird things happen when P --> 0 TODO FIXME
            
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_trials))';
            multi.durations{1} = zeros(size(data.contextRole(which_trials)));
            
            multi.pmod(1).name{1} = 'surprise';
            multi.pmod(1).param{1} = surprise(6:end)';
            multi.pmod(1).poly{1} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_trials))';
            multi.durations{2} = zeros(size(data.contextRole(which_trials)));
            
            
        % M4 & M1 value pmod @ trial onset, before updated
        % + M4 & M1 update pmods @ feedback (outcome) onset, after update
        % result: after clust FWE, M4_value overlaps with posterior hippo a bit but it's mostly
        %         temporal / insula + S1
        %         M4_update is all over the place
        %
        case 104
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M4_value';
            multi.pmod(1).param{1} = valuess(:,4)'; % expected outcome predicted by M4 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M4 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M4_update';
            multi.pmod(2).param{1} = new_valuess(:,4)' - valuess(:,4)'; % PE = new expected outcome - old expected outcome by M4 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = new_valuess(:,1)' - valuess(:,1)'; % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order
            
        % M4 & M1 value pmod @ trial onset, before updated
        % + M4 & M1 update pmods @ feedback (outcome) onset, after update
        % same as 104 but with absolute updates
        % result: nothing really; M4_value is kind of temporal again;
        %         M4_update - M1_update is negative in posterior hippo
        %
        case 105
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M4_value';
            multi.pmod(1).param{1} = valuess(:,4)'; % expected outcome predicted by M4 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M4 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M4_update';
            multi.pmod(2).param{1} = abs(new_valuess(:,4)' - valuess(:,4)'); % PE = new expected outcome - old expected outcome by M4 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = abs(new_valuess(:,1)' - valuess(:,1)'); % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order

        % M2 value pmod @ trial onset, before updated
        % + M2 update pmod @ feedback (outcome) onset, after update
        % similar to 105 but with absolute changes
        % result: M2_value with posterior hippo; kinda like M4_value in 104
        %         but smaller; does NOT survive cluster FWE
        %         M2_update overlaps with striatum but also with everything
        %         else
        %
        case 106
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M2_value';
            multi.pmod(1).param{1} = valuess(:,2)'; % expected outcome predicted by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % M2 update ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M2_update';
            multi.pmod(2).param{1} = abs(new_valuess(:,2)' - valuess(:,2)'); % PE = new expected outcome - old expected outcome by M2 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));            
            
        % value pmod @ trial onset, before updated
        % + update pmod @ feedback (outcome) onset, after update
        % similar to 106 but with total value
        % result: value in posterior hippo; not after FWE correction
        %         BUT ...very suspicious of the negative activations there
        %         at feedback time
        %
        case 107
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'value';
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % update ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'update';
            multi.pmod(2).param{1} = abs(new_values' - values'); % PE = new expected outcome - old expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));            

        % value pmod @ trial onset, before updated
        % + abs update pmod @ feedback (outcome) onset, after update
        % same as 107 but total update
        % result: similar to 107
        %
        case 108
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'value';
            multi.pmod(1).param{1} = values'; % expected outcome for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order
            
            % update ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'update';
            multi.pmod(2).param{1} = new_values' - values'; % PE = new expected outcome - old expected outcome for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));            

            
        % M2 & M1 value pmod @ trial onset, before updated
        % + M2 & M1 update pmods @ feedback (outcome) onset, after update
        % same idea as 104 except M2 instead of M4
        % result: nothing
        %
        case 109
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M2_value';
            multi.pmod(1).param{1} = valuess(:,2)'; % expected outcome predicted by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M2 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M2_update';
            multi.pmod(2).param{1} = new_valuess(:,2)' - valuess(:,2)'; % PE = new expected outcome - old expected outcome by M2 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = new_valuess(:,1)' - valuess(:,1)'; % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order

        % M3 & M1 value pmod @ trial onset, before updated
        % + M3 & M1 update pmods @ feedback (outcome) onset, after update
        % same idea as 104 except M3 instead of M4
        % result: M3_value -> bilateral posterior hippocampus activation;
        %                     doesn't survive cluster FWE
        %
        case 110
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M3_value';
            multi.pmod(1).param{1} = valuess(:,3)'; % expected outcome predicted by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M3 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M3_update';
            multi.pmod(2).param{1} = new_valuess(:,3)' - valuess(:,3)'; % PE = new expected outcome - old expected outcome by M3 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = new_valuess(:,1)' - valuess(:,1)'; % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order
            

        % M2 & M1 value pmod @ trial onset, before updated
        % + M2 & M1 update pmods @ feedback (outcome) onset, after update
        % same as 109 except with const reaction regressor
        % result: nothing
        %
        case 111
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M2_value';
            multi.pmod(1).param{1} = valuess(:,2)'; % expected outcome predicted by M2 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M2 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M2_update';
            multi.pmod(2).param{1} = new_valuess(:,2)' - valuess(:,2)'; % PE = new expected outcome - old expected outcome by M2 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = new_valuess(:,1)' - valuess(:,1)'; % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order

            % const @ RT (trials 1..20)
            % 
            multi.names{3} = 'RT';
            multi.onsets{3} = -1 + cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));            
            
        % M3 & M1 value pmod @ trial onset, before updated
        % + M3 & M1 update pmods @ feedback (outcome) onset, after update
        % same as 110 except with const reaction regressor
        % result: same as 110
        %
        case 112
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualChoiceOnset(which_train))';
            multi.durations{1} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(1).name{1} = 'M3_value';
            multi.pmod(1).param{1} = valuess(:,3)'; % expected outcome predicted by M3 for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order

            multi.pmod(1).name{2} = 'M1_value';
            multi.pmod(1).param{2} = valuess(:,1)'; % expected outcome predicted by M1 for trials 1..20
            multi.pmod(1).poly{2} = 1; % first order
            
            % M1 & M3 updates ~= prediction errors @ feedback
            %
            multi.names{2} = 'outcome';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{2} = ones(size(data.contextRole(which_train))); % 1 s durations
            
            multi.pmod(2).name{1} = 'M3_update';
            multi.pmod(2).param{1} = new_valuess(:,3)' - valuess(:,3)'; % PE = new expected outcome - old expected outcome by M3 for trials 1..20
            multi.pmod(2).poly{1} = 1; % first order
            
            multi.pmod(2).name{2} = 'M1_update';
            multi.pmod(2).param{2} = new_valuess(:,1)' - valuess(:,1)'; % PE = new expected outcome - old expected outcome by M1 for trials 1..20
            multi.pmod(2).poly{2} = 1; % first order
            
            % const @ RT (trials 1..20)
            % 
            multi.names{3} = 'RT';
            multi.onsets{3} = -1 + cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));            
            
            
            
        % Prediction error using subject choices @ feedback
        % result: 'actual' -- weird left cuneus after cluster FWE
        %                    for P = 0.005, also mPFC
        %
        case 113
            % Prediction error @ feedback
            %
            multi.names{1} = 'prediction_error';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            pressed_sick = strcmp(data.response.keys(which_train), 'left');
            PE = outcomes' - pressed_sick'; % outcome - subject predicted outcome for trials 1..20
            if std(PE) > 1e-6 % only include if variable
                multi.pmod(1).name{1} = 'actual';
                multi.pmod(1).param{1} = PE;
                multi.pmod(1).poly{1} = 1; % first order  
            end
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
                        
        % Prediction error using subject choices @ feedback
        % same as 113 + const regressor @ RT
        % WRONG: doesn't work, also events too close together
        %
        case 114
            % Prediction error @ feedback
            %
            multi.names{1} = 'prediction_error';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            pressed_sick = strcmp(data.response.keys(which_train), 'left');
            PE = outcomes' - pressed_sick'; % outcome - subject predicted outcome for trials 1..20
            if std(PE) > 1e-6 % only include if variable
                multi.pmod(1).name{1} = 'actual';
                multi.pmod(1).param{1} = PE;
                multi.pmod(1).poly{1} = 1; % first order  
            end
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % const @ RT (trials 1..20) TODO some trials don't have
            % responses
            %
            multi.names{3} = 'keypress';
            multi.onsets{3} = -1 + cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));

            
        % regressors 1s around trial onset time
        % to see how things change around it
        % results:
        %    before_trial_onset: low visual; neg motor & medial
        %         frontal
        %    trial_onset: lots visual, left motor (right hand)!; neg
        %         medial frontal
        %    after_trial_onset: visual, neg med frontal
        %    trial_onset - before_trial_onset: left motor (right
        %         hand), med cinculate cortex; NO visual?
        %    after_trial_onset - before_trial_onset: lots visual!, left
        %         motor, ACC
        %    after_trial_onset - trial_onset: lots mPFC!, neg visual!
        % 
        case 115
            % const @ trial onset - 1 (trials 1..24)
            % 
            multi.names{1} = 'before_trial_onset';
            multi.onsets{1} = -1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));

            % const @ trial onset (trials 1..24)
            % 
            multi.names{2} = 'at_trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset + 1 (trials 1..24)
            % 
            multi.names{3} = 'after_trial_onset';
            multi.onsets{3} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows)));
            
        % regressors 0.5s around trial onset time
        % to see how things change around it
        % similar to 115
        % results: poor version of 115; ignore these (too close)
        %
        case 116
            % const @ trial onset - 0.5 (trials 1..24)
            % 
            multi.names{1} = 'before_trial_onset';
            multi.onsets{1} = -0.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));

            % const @ trial onset (trials 1..24)
            % 
            multi.names{2} = 'at_trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset + 0.5 (trials 1..24)
            % 
            multi.names{3} = 'after_trial_onset';
            multi.onsets{3} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows)));

        % regressors before trial onset time
        % to see how things change around it
        % similar to 116
        % results: 
        %    1.5_before_trial_onset: pos visual & left motor?; neg mPFC?
        %    1_before_trial_onset: neg visual & left motor; pos mPFC ?
        %                   (opposite of above)
        %    0.5_before_trial_onset: same as first one
        %    0.5_before_trial_onset - 1.5_before_trial_onset: lots of very
        %              high activity in visual & parietal
        %    1_before_trial_onset - 1.5_before_trial_onset: neg visual &
        %              left motor, pos mPFC
        %    0.5_before_trial_onset - 1_before_trial_onset: negative of
        %              above
        %
        case 117
            % const @ trial onset - 1.5 (trials 1..24)
            % 
            multi.names{1} = '1.5_before_trial_onset';
            multi.onsets{1} = -1.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));

            % const @ trial onset - 1 (trials 1..24)
            % 
            multi.names{2} = '1_before_trial_onset';
            multi.onsets{2} = -1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset - 0.5 (trials 1..24)
            % 
            multi.names{3} = '0.5_before_trial_onset';
            multi.onsets{3} = -0.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows)));

        % regressors after trial onset time
        % to see how things change around it
        % similar to 116
        % results:
        %   1_after_trial_onset - 0.5_after_trial_onset: neg left motor;
        %     low pos visual
        %   1_after_trial_onset - trial_onset: pos mPFC; neg left motor &
        %   visual
        %   0.5_after_trial_onset - trial_onset: same
        %
        case 118
            % const @ trial onset (trials 1..24)
            % 
            multi.names{1} = '0_trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));

            % const @ trial onset + 0.5 (trials 1..24)
            % 
            multi.names{2} = '0.5_after_trial_onset';
            multi.onsets{2} = 0.5 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
            % const @ trial onset + 1 (trials 1..24)
            % 
            multi.names{3} = '1_after_trial_onset';
            multi.onsets{3} = 1 + cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{3} = zeros(size(data.contextRole(which_rows)));
            
        % regressors 0.5s around feedback time
        % to see how things change around it
        % similar to 115
        % results:
        %   feedback_onset - before_feedback_onset: pos mPFC; neg left
        %      motor; low visual
        %   after_feedback_onset - before_feedback_onset: massive neg
        %      visual, left motor & cingulate; pos mPFC
        %   afterfeedback_onset - feedback_onset: neg visual & motor; pos
        %      mPFC
        %
        case 119
            % const @ feedback onset - 0.5 (trials 1..20)
            % 
            multi.names{1} = 'before_feedback_onset';
            multi.onsets{1} = -0.5 + cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            % const @ feedback onset (trials 1..20)
            % 
            multi.names{2} = 'at_feedback_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % const @ feedback onset + 0.5 (trials 1..20)
            % 
            multi.names{3} = 'after_feedback_onset';
            multi.onsets{3} = 0.5 + cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{3} = zeros(size(data.contextRole(which_train)));
            

        % continuous regressors during iti and trials
        % results:
        %   trial = positive activations everywhere; visual in particular
        %   iti = negative everywhere
        %   trial - iti: positive everywhere (!!YES!!!)
        %
        case 120
            trial_onsets = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            iti_onsets = cellfun(@str2num, data.actualItiOnset(which_rows))';
            iti_offsets = cellfun(@str2num, data.actualItiOffset(which_rows))';
            
            % const during trials (trials 1..24)
            % 
            multi.names{1} = 'trial';
            multi.onsets{1} = trial_onsets;
            multi.durations{1} = iti_onsets - trial_onsets;

            % const during ITIs (trials 1..24)
            % 
            multi.names{2} = 'iti';
            multi.onsets{2} = iti_onsets;
            multi.durations{2} = iti_offsets - iti_onsets;
       
        % main effect @ trial onset, test trials only
        % result: nothing
        %
        case 121
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_test))';
            multi.durations{1} = zeros(size(data.contextRole(which_test)));
            
            
        % correct vs. wrong pmod @ outcome
        % similar to 3 but done right
        % result: 'wrong': anterior insula, mPFC (after FWE)
        %
        case 122
            which_error = which_train & ~data.response.corr;
            
            % const @ feedback onset (trials 1..20)
            %
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
            
        % Bayesian surprise = Kullback?Leibler divergence @ feedback
        % (outcome) onset
        % + error regressor on wrong trials, to account for error signal
        % result: THIS IS IT -- R angular gyrus, lateral OFC, R lateral PFC
        %          (after FWE)
        % #KEEP in dropbox
        %
        case 123
            which_error = which_train & ~data.response.corr;
            
            priors = which_structures / sum(which_structures);
            Q = [priors; P(1:end-1,:)];
            logs = log2(P) - log2(Q); 
            logs(isnan(logs)) = 0; % lim_{x->0} x log(x) = 0
            surprise = sum(P .* logs, 2);
            surprise(isnan(surprise)) = 0; % weird things happen when P --> 0, e.g. we get -Infs

            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'surprise';
            multi.pmod(1).param{1} = surprise';
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
            
        % main effect BUT with the most likely structure (not the actual
        % condition) based on subject's choices
        % result: nothing
        %
        case 124
            % sanity check
            %{ 
            x1c3_choice = data.response.keys(which_test & data.cueId == 0 & data.contextId == 2);
            x3c1_choice = data.response.keys(which_test & data.cueId == 2 & data.contextId == 0);
            predict = @(V_n) softmax(V_n, inv_softmax_temp);
            if ~strcmp(x1c3_choice, 'None') && ~strcmp(x3c1_choice, 'None')
                x1c3_choice = strcmp(x1c3_choice, 'left');
                x3c1_choice = strcmp(x3c1_choice, 'left');
                
                x1c3_idx = data.trialId(which_test & data.cueId == 0 & data.contextId == 2);
                x3c1_idx = data.trialId(which_test & data.cueId == 2 & data.contextId == 0);
                
                x1c3_values = test_valuess(x1c3_idx, :);
                x3c1_values = test_valuess(x3c1_idx, :);
               
                x1c3_probs = predict(x1c3_values);
                x3c1_probs = predict(x3c1_values);
                
                x1c3_likelihoods = binopdf(x1c3_choice, 1, x1c3_probs);
                x3c1_likelihoods = binopdf(x3c1_choice, 1, x3c1_probs);
            end
            %}
            
            % for sanity check with above
            %which_trials = which_test & data.cueId ~= data.contextId & ~strcmp(data.response.keys, 'None');
            
            which_trials = which_test & ~strcmp(data.response.keys, 'None');
            if sum(which_trials) > 0
                subj_choices = strcmp(data.response.keys(which_trials), 'left');
                trial_idxs = data.trialId(which_trials);

                trial_values = test_valuess(trial_idxs, :);
                predict = @(V_n) softmax(V_n, inv_softmax_temp);
                trial_probs = predict(trial_values);

                %trial_likelihoods = arrayfun(@(i) binopdf(subj_choices(i), 1, trial_probs(i,:)), 1:numel(subj_choices), 'UniformOutput', false);
                %trial_likelihoods = reshape(cell2mat(trial_likelihoods'), size(trial_probs)); % careful -- must transpose first
                % FUCK ncf...
                liks = nan(size(trial_probs));
                for i = 1:size(trial_probs, 1) % for each test trial
                    for j = 1:size(trial_probs, 2) % for each causal structure
                        liks(i,j) = binopdf(subj_choices(i), 1, trial_probs(i, j));
                    end
                end
                %assert(immse(liks, trial_likelihoods) < 1e-10);

                % take average log likelihood !!!
                % this is to account for missing trials where subject timed
                % out. Note that taking the sum would be wrong -- imagine
                % if subject responded on only 1 trial
                test_log_lik = mean(log(liks), 1);
                test_log_lik = test_log_lik(:, 1:3); % exclude M4
                [~, most_likely_M] = max(test_log_lik);

                % see if we match the actual condition, just FYI
                %
                M = -1;
                if strcmp(condition, 'irrelevant')
                    M = 1;
                elseif strcmp(condition, 'modulatory')
                    M = 2;
                else
                    assert(strcmp(condition, 'additive'));
                    M = 3;
                end
                fprintf('   M = %d, most likely M = %d\n', M, most_likely_M);
                if M ~= most_likely_M
                    disp('                DIFFERENT!');
                end

                % most likely causal structure @ feedback
                %
                switch most_likely_M
                    case 1
                        multi.names{1} = 'irrelevant';
                    case 2
                        multi.names{1} = 'modulatory';
                    case 3
                        multi.names{1} = 'additive';
                end
                multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
                multi.durations{1} = zeros(size(data.contextRole(which_train)));

                % const @ trial onset (trials 1..20)
                % 
                multi.names{2} = 'trial_onset';
                multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
                multi.durations{2} = zeros(size(data.contextRole(which_train)));
            end
            
        % Cause-specific posterior @ feedback time
        %
        case 125
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            multi.pmod(1).name{1} = 'posterior';
            multi.pmod(1).param{1} = P(:,M)'; % cause-specific posterior for trials 1..20
            multi.pmod(1).poly{1} = 1; % first order  

            % const @ trial onset
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % Main effect boxcar for whole trial @ trial onset
        % result: nothing; the 'additive - irrelevant' hippocampus thing
        % disappears
        %
        case 126
            % context role @ feedback/outcome onset
            % 
            trial_onsets = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            iti_onsets = cellfun(@str2num, data.actualItiOnset(which_train))';
            multi.names{1} = condition;
            multi.onsets{1} = trial_onsets;
            multi.durations{1} = iti_onsets - trial_onsets;
            
            
        %
        % These are with the random effects parameters fit for the fMRI
        % behavioral data
        %
            
        % Bayesian surprise = Kullback?Leibler divergence @ feedback
        % (outcome) onset
        % + error regressor on wrong trials, to account for error signal
        % same as 123
        %
        case 127     
            surprise = simulated.surprise(which_train, :);
    
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'surprise';
            multi.pmod(1).param{1} = surprise';
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            which_error = which_train & ~data.response.corr;
            
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
         
        % Squared PE = outcome - expected outcome @ feedback time
        %
        case 128
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'squared_PE';
            multi.pmod(1).param{1} = (outcomes' - values').^2;
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            
        % Squared PE 2 = new expected outcome - expected outcome @ feedback time
        %            
        case 129
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'squared_PE_2';
            multi.pmod(1).param{1} = (new_values' - values').^2;
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M1 posterior
        %            
        case 130
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_posterior';
            multi.pmod(1).param{1} = P(:, 1)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));


        % M2 posterior
        %            
        case 131
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_posterior';
            multi.pmod(1).param{1} = P(:, 2)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            

        % M3 posterior
        %            
        case 132
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_posterior';
            multi.pmod(1).param{1} = P(:, 3)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            
        % M1 reliability or precision = 1/lambda
        %            
        case 133
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M1_reliability';
            multi.pmod(1).param{1} = 1./lambdas(:, 1)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % M2 reliability or precision = 1/lambda
        %            
        case 134
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M2_reliability';
            multi.pmod(1).param{1} = 1./lambdas(:, 2)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            

        % M3 reliability or precision = 1/lambda
        %            
        case 135
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'M3_reliability';
            multi.pmod(1).param{1} = 1./lambdas(:, 3)';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % main effect @ trial onset
        % compare 141 (feedback) with 136 (trial onset) with 138 (null)
        %
        case 136
            multi.names{1} = 'feedback_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.names{2} = condition;
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            

        % exactly the same as 123 -- sanity check the new code
        %
        case 137
            which_error = which_train & ~data.response.corr;
            
            priors = which_structures / sum(which_structures);
            Q = [priors; P(1:end-1,:)];
            logs = log2(P) - log2(Q); 
            logs(isnan(logs)) = 0; % lim_{x->0} x log(x) = 0
            surprise = sum(P .* logs, 2);
            surprise(isnan(surprise)) = 0; % weird things happen when P --> 0, e.g. we get -Infs

            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'surprise';
            multi.pmod(1).param{1} = surprise';
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
            
        % null GLM @ trial onset
        % compare 141 (feedback) with 136 (trial onset) with 138 (null)
        % result: 'feedback_onset - trial_onset' NICE BILATERAL
        %          HIPPOCAMPUS!!!
        %
        case 138
            multi.names{1} = 'feedback_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            

        % null GLM @ trial onset (TRAINING trials only)
        % compare 1 (feedback) with 140 (trial onset) with 139 (null)
        %
        case 139
            multi.names{1} = 'feedback_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            
        % main effect @ trial onset (TRAINING trials only)
        % compare 1 (feedback) with 140 (trial onset) with 139 (null)
        %
        case 140
            multi.names{1} = 'feedback_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.names{2} = condition;
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
        % main effect @ feedback onset
        % compare 141 (feedback) with 136 (trial onset) with 138 (null)
        %
        case 141
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{2} = zeros(size(data.contextRole(which_rows)));
            
            
        % correct vs wrong @ feedback
        %
        case 142
            which_error = which_train & ~data.response.corr;
            which_correct = which_train & ~which_error;
            assert(sum(which_error) + sum(which_correct) == 20);
            
            multi.names{1} = 'trial_onset';
            multi.onsets{1} = cellfun(@str2num, data.actualChoiceOnset(which_rows))';
            multi.durations{1} = zeros(size(data.contextRole(which_rows)));

            multi.names{2} = 'correct';
            multi.onsets{2} = cellfun(@str2num,data.actualFeedbackOnset(which_correct))';
            multi.durations{2} = zeros(size(data.contextRole(which_correct)));
            
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
        % For each trial, have a trial_onset regressor and a feedback regressor (training trials only)
        % it's for the classifier
        % #KEEP in dropbox
        %
        case 143
            idx = 0;

            trial_onsets = data.actualChoiceOnset(which_rows);
            for t=1:metadata.trialsPerRun
                idx = idx + 1;
                multi.names{idx} = ['trial_onset_', num2str(t)];
                multi.onsets{idx} = [str2double(trial_onsets(t))];
                multi.durations{idx} = [0];
            end
            
            feedback_onsets = data.actualFeedbackOnset(which_train);
            for t=1:metadata.trainingTrialsPerRun
                idx = idx + 1;
                multi.names{idx} = ['feedback_onset_', num2str(t)];
                multi.onsets{idx} = [str2double(feedback_onsets(t))];
                multi.durations{idx} = [0];
            end

        % KL for structures vs. KL for weights
        % where do we do parameter estimation vs. causal inference
        % result: KL_structures same as before, KL_weights all over the place (but see 145)
        %
        case 144
            which_error = which_train & ~data.response.corr;
            
            KL_structures = simulated.surprise(which_train);
            KL_weights = simulated.KL_weights(which_train);

            % put the regressors
            %
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'KL_weights';
            multi.pmod(1).param{1} = KL_weights';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'KL_structures';
            multi.pmod(1).param{2} = KL_structures';
            multi.pmod(1).poly{2} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
       
        % KL for sturcts vs. KL for weights
        % same as 144 but with error regressor
        % result: KL_structures same as before, KL_weights - wrong in right S1/M1 only
        %
        case 145
            which_error = which_train & ~data.response.corr;
            
            KL_structures = simulated.surprise(which_train);
            KL_weights = simulated.KL_weights(which_train);

            % put the regressors
            %
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'KL_weights';
            multi.pmod(1).param{1} = KL_weights';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'KL_structures';
            multi.pmod(1).param{2} = KL_structures';
            multi.pmod(1).poly{2} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
        % KL for structures vs. KL for weights
        % same as 144 but swapped structure and weights
        % result: KL_structures same as before, KL_weights all over the place (but see 147)
        %
        case 146
            which_error = which_train & ~data.response.corr;
            
            KL_structures = simulated.surprise(which_train);
            KL_weights = simulated.KL_weights(which_train);

            % put the regressors
            %
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights';
            multi.pmod(1).poly{2} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
       
        % KL for sturcts vs. KL for weights
        % same as 145 but swapped structure and weights
        % result: KL_structures same as before, KL_weights - wrong in right S1/M1 only
        %
        case 147
            which_error = which_train & ~data.response.corr;
            
            KL_structures = simulated.surprise(which_train);
            KL_weights = simulated.KL_weights(which_train);

            % put the regressors
            %
            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num,data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order        

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights';
            multi.pmod(1).poly{2} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));
            
            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
        
        % KL for weights for strucutre corresponding to condition
        % THIS IS IT -- going in paper
        % #KEEP in dropbox
        %
        case 148
            which_error = which_train & ~data.response.corr;

            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            ww_prior = simulated.ww_before{M}(which_train, :);
            Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
            ww_posterior = simulated.ww_after{M}(which_train, :);
            Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
            KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
            assert(size(KL_weights, 1) == sum(which_train));
            
            multi.pmod(1).name{1} = 'KL_weights';
            multi.pmod(1).param{1} = KL_weights';
            multi.pmod(1).poly{1} = 1; % first order        

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            

        % KL for weights for all conditions, not orthogonalized
        %
        case 149
            which_error = which_train & ~data.response.corr;

            multi.names{1} = 'feedback';
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            multi.orth{1} = 0; % do NOT orthogonalize the feedback event

            for M = 1:3
                ww_prior = simulated.ww_before{M}(which_train, :);
                Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
                ww_posterior = simulated.ww_after{M}(which_train, :);
                Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
                KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
                assert(size(KL_weights, 1) == sum(which_train));
                
                multi.pmod(1).name{M} = ['KL_weights_M', num2str(M)];
                multi.pmod(1).param{M} = KL_weights';
                multi.pmod(1).poly{M} = 1; % first order        
            end

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
            
        % KL for weights for strucutre corresponding to condition
        % + KL for posterior
        % similar to 148
        %
        case 150
            which_error = which_train & ~data.response.corr;

            KL_structures = simulated.surprise(which_train);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            ww_prior = simulated.ww_before{M}(which_train, :);
            Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
            ww_posterior = simulated.ww_after{M}(which_train, :);
            Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
            KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
            assert(size(KL_weights, 1) == sum(which_train));
            
            multi.pmod(1).name{1} = 'KL_weights';
            multi.pmod(1).param{1} = KL_weights';
            multi.pmod(1).poly{1} = 1; % first order        
            

            multi.pmod(1).name{2} = 'KL_structures';
            multi.pmod(1).param{2} = KL_structures';
            multi.pmod(1).poly{2} = 1; % first order                    

            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
        % Same as 150 but swapped KL_weights and KL_posterior
        %
        case 151
            which_error = which_train & ~data.response.corr;

            KL_structures = simulated.surprise(which_train);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            ww_prior = simulated.ww_before{M}(which_train, :);
            Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
            ww_posterior = simulated.ww_after{M}(which_train, :);
            Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
            KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
            assert(size(KL_weights, 1) == sum(which_train));
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order                    

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights';
            multi.pmod(1).poly{2} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
            
        % Same as 151 but averaged KL_weights across the 3 structures,
        % weighted by the structure prior
        %
        case 152
            which_error = which_train & ~data.response.corr;

            KL_structures = simulated.surprise(which_train);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            KL_weights = [];
            for M = 1:3
                ww_prior = simulated.ww_before{M}(which_train, :);
                Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
                ww_posterior = simulated.ww_after{M}(which_train, :);
                Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
                KL_weights(:, M) = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
                assert(size(KL_weights, 1) == sum(which_train));
            end
            structure_priors = simulated.Q(which_train, 1:3);
            KL_weights_avg = sum(KL_weights .* structure_priors, 2);
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order                    

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights_avg';
            multi.pmod(1).poly{2} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
            
        % Same as 151 but with orthogonalization off
        %
        case 153
            which_error = which_train & ~data.response.corr;

            KL_structures = simulated.surprise(which_train);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            multi.orth{1} = 0; % do NOT orthogonalize them!
                        
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            ww_prior = simulated.ww_before{M}(which_train, :);
            Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
            ww_posterior = simulated.ww_after{M}(which_train, :);
            Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
            KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
            assert(size(KL_weights, 1) == sum(which_train));
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order                    

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights';
            multi.pmod(1).poly{2} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
            
            
        % Same as 153 but with context changes as a regressor
        %
        case 154
            which_error = which_train & ~data.response.corr;

            KL_structures = simulated.surprise(which_train);
            
            context_changed = data.contextId ~= circshift(data.contextId, 1);
            context_changed(data.trialId == 1 & data.isTrain) = 0; % TODO 0 or 1?
            context_changed = context_changed(which_train);
            
            % context role @ feedback/outcome onset
            % 
            multi.names{1} = condition;
            multi.onsets{1} = cellfun(@str2num, data.actualFeedbackOnset(which_train))';
            multi.durations{1} = zeros(size(data.contextRole(which_train)));

            multi.orth{1} = 0; % do NOT orthogonalize them!
                        
            M = -1;
            if strcmp(condition, 'irrelevant')
                M = 1;
            elseif strcmp(condition, 'modulatory')
                M = 2;
            else
                assert(strcmp(condition, 'additive'));
                M = 3;
            end
            ww_prior = simulated.ww_before{M}(which_train, :);
            Sigma_prior = simulated.Sigma_before{M}(:,:,which_train);
            ww_posterior = simulated.ww_after{M}(which_train, :);
            Sigma_posterior = simulated.Sigma_after{M}(:,:,which_train);
            KL_weights = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);
            assert(size(KL_weights, 1) == sum(which_train));
            
            multi.pmod(1).name{1} = 'KL_structures';
            multi.pmod(1).param{1} = KL_structures';
            multi.pmod(1).poly{1} = 1; % first order                    

            multi.pmod(1).name{2} = 'KL_weights';
            multi.pmod(1).param{2} = KL_weights';
            multi.pmod(1).poly{2} = 1; % first order        
            
            % const @ trial onset (trials 1..20)
            % 
            multi.names{2} = 'trial_onset';
            multi.onsets{2} = cellfun(@str2num, data.actualChoiceOnset(which_train))';
            multi.durations{2} = zeros(size(data.contextRole(which_train)));

            % correct vs. wrong (0/1) @ feedback / outcome onset (WRONG trials 1..20)
            % 
            if sum(which_error) > 0
                multi.names{3} = 'wrong';
                multi.onsets{3} = cellfun(@str2num,data.actualFeedbackOnset(which_error))';
                multi.durations{3} = zeros(size(data.contextRole(which_error)));
            end
                        
            % did context change on this trial?
            %
            multi.names{4} = 'context_changed';
            multi.onsets{4} = cellfun(@str2num,data.actualFeedbackOnset(context_changed))';
            multi.durations{4} = zeros(size(data.contextRole(context_changed)));
            

            
        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');
            
    
    end % end of switch statement

   if save_output
       save('context_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end 
