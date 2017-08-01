function [test_log_liks, test_RTs] = get_test_behavior()
% Get the subject test trial behavior
%
% INPUT: n/a
%
% OUTPUT:
% test_log_liks = N x R matrix of average test choice log likelihoods according to
%                 the model, where N = # of subjects, R = # runs 
% test_RTs = N x R matrix of average test choice RTs
%


%% Load behavior
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
goodSubjects = getGoodSubjects();
which_rows = data.which_rows;

%% Simulate behavior
%
prior_variance = 0.1249;
inv_softmax_temp = 2.0064;
params = [prior_variance inv_softmax_temp];
which_structures = logical([1 1 1 0]);
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows, false);

%% Get test choice log likelihood and RTs
%
test_log_liks = nan(metadata.N, metadata.runsPerSubject);
test_RTs = nan(metadata.N, metadata.runsPerSubject);

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

        % Get the test trial likelihood
        %
        run_test_trials = run_trials & ~data.isTrain & ~data.timeout;
        
        X_fixed = data.chose_sick(run_test_trials); % actual subject choice on each trial
        P = simulated.pred(run_test_trials); % probability of subject choice on each trial
        assert(numel(X_fixed) <= metadata.testTrialsPerRun);
        
        liks = binopdf(X_fixed, 1, P);
        assert(numel(liks) == numel(X_fixed));
        % average to account for timeouts
        % this is important -- summing them would not be okay e.g. if
        % subject responded on only 1 trial
        %
        avg_loglik = mean(log(liks));
        %avg_loglik = prod(liks); DON'T DO IT
        
        test_log_liks(subj_idx, run) = avg_loglik; 
        
        % Get the test trial RTs
        %
        test_RTs(subj_idx, run) = mean(data.response.rt(run_test_trials));
    end
end





