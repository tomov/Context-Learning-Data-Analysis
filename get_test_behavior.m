function [test_log_liks, test_RTs, test_log_liks_all] = get_test_behavior(params, which_structures)
% Get the subject test trial behavior
%
% INPUT: 
% [params, which_structures] = model_default_params();
%
% OUTPUT:
% test_log_liks = N x R matrix of average test choice log likelihoods according to
%                 the model, where N = # of subjects, R = # runs 
% test_RTs = N x R matrix of average test choice RTs
% test_log_liks_all = N x R x 3 matrix of average test choice log likelihoods for each structure separately
%


%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Simulate behavior
%
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows, false);

%predict = @(V_n) softmax(V_n, inv_softmax_temp);

%% Get test choice log likelihood and RTs
%
test_log_liks = nan(metadata.N, metadata.runsPerSubject);
test_RTs = nan(metadata.N, metadata.runsPerSubject);
test_log_liks_all = nan(metadata.N, metadata.runsPerSubject, sum(which_structures));

% Iterate over subjects
%
subj_idx = 0; % 1..20
for subj = goodSubjects % 1..25
    subject = metadata.allSubjects(subj); % 'con001' ... 'con025'
    subj_trials = data.which_rows & strcmp(data.participant, subject);
    subj_idx = subj_idx + 1;
    fprintf('subj %d (idx %d)\n', subj, subj_idx);

    % Iterate over runs
    %
    for run = 1:metadata.runsPerSubject
        fprintf(' RUN %d\n', run);
        run_trials = subj_trials & data.runId == run;
        condition = data.contextRole(run_trials);
        condition = condition{1};
		assert(sum(run_trials) == metadata.trialsPerRun);

        % Get the overall test trial likelihood
        %
        run_test_trials = run_trials & ~data.isTrain & ~data.timeout;
        
        X_fixed = data.chose_sick(run_test_trials); % actual subject choice on each trial
        p = simulated.pred(run_test_trials); % probability of subject choice on each trial
        assert(numel(X_fixed) <= metadata.testTrialsPerRun);
        
        liks = binopdf(X_fixed, 1, p);
        assert(numel(liks) == numel(X_fixed));
        % average to account for timeouts
        % this is important -- summing them would not be okay e.g. if
        % subject responded on only 1 trial
        %
        avg_loglik = mean(log(liks));
        %avg_loglik = prod(liks); DON'T DO IT
        
        test_log_liks(subj_idx, run) = avg_loglik; 

        %{
        % Get the test trial likelihood according to each structure
        %
        v = simulated.valuess(run_test_trials, which_structures); % value on each trial, for each structure
        p = predict(v); % probability of subject choice on each trial, for each structure
        for s = 1:size(p,2) % for each structure, get the average log likelihood
            liks = binopdf(X_fixed, 1, p(:,s));
            assert(numel(liks) == numel(X_fixed));

            avg_loglik = mean(log(liks));
            test_log_liks_all(subj_idx, run, s) = avg_loglik;
        end
        %}
        
        % Get the test trial RTs
        %
        test_RTs(subj_idx, run) = mean(data.response.rt(run_test_trials));
    end
end





