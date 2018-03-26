% Idea: attention (multivariate) <-> condition:
% 1) for a given run,
% - in the irrelevant context condition, you only need to attend to the cue 
% => regions sensitive to that would look similar in trials with the same cue, event when context is different 
% - conversely for the additive (irrelevant cue) condition
% - for modulatory, both cue & context matter 
% so use a region's similarity with a cue-only RDM, context-only RDM, or both RDM as a (neural) proxy for whether that region "thinks"/"prefers" structures M1, M1', or M2 
% use spheres around the peak voxels
% 2) separately, compute P(M1 | test choices) = P(test choices | M1) P(M1) / const = test lik under M1 / const 
% do it for all structures
% this tells us which structure the subject's test choices are most consistent with 
% 3) collect those from all runs for a given subject => two 27-element vectors (9 runs * 3 numbers per run):
% - one with the spearman rho's
% - one with the posteriors over structures based on test choices
% and correlate them
%
% kinda borrows from rdms_behavior.m and correlate_neural_and_behavior.m
%

utils;

dirname = 'rdms';

[params, which_structures] = model_params('results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat'); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! important which one! !!!!!!!!
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

which_trials = data.which_rows & data.isTrain; % Look at training trials only

%% get the neural rdms
%
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
%r = 1.814;
r = 2.6667;  
%r = 6.6667; 
direct = '+';

use_nosmooth = false;

Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 171, 'KL_structures', p, direct, alpha, Dis, Num, r, {'trial_onset', 'feedback_onset'}, use_nosmooth); %  <------ nothing :( nosmooth: p = 0.009 for angular_R at trial_onset; p = 0.03 for Parietal_Inf at feedback onset
%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, 'rdms/M1M2M1_4mm/searchlight_tmap_posterior_feedback_onset.nii', 0, 'light', p, direct, alpha, Dis, Num, r, {'trial_onset', 'feedback_onset'}, use_nosmooth);  % <-- nothing :( nosmooth: some at p = 0.01, 0.03, 0.04 but won't survive...
%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, 'rdms/M1M2M1_4mm/searchlight_tmap_prior_trial_onset.nii', 0, 'light', p, direct, alpha, Dis, Num, r, {'trial_onset', 'feedback_onset'}, use_nosmooth); % <--  nothing
%Neural = rdms_get_spheres_from_contrast(data, metadata, which_trials, context_expt(), 170, 'KL_clusters', p, direct, alpha, Dis, Num, r, {'trial_onset', 'feedback_onset'}, use_nosmooth); %  <---  nothing



%% find how similar representations are to posterior (model_idx = 1) in different ROIs per-run
%
M1_rhos = nan(metadata.N, metadata.runsPerSubject, numel(Neural));
M2_rhos = nan(metadata.N, metadata.runsPerSubject, numel(Neural));
M4_rhos = nan(metadata.N, metadata.runsPerSubject, numel(Neural)); % technically M1'

goodSubjects = getGoodSubjects();
subjs = metadata.allSubjects(goodSubjects);

% t1, t2 are trial indices corresponding to the pairs of trials in each
% cell of the RDM (for a single subject)
% used to quickly generate RDM masks and extract sub-RDMs.
% Assumes the same trials were used from all subjects (so it takes subj 1
% for convenience)
%
which_trials_per_subj = which_trials & strcmp(data.participant, subjs{1});
[t1, t2] = meshgrid(find(which_trials_per_subj), find(which_trials_per_subj));


[cueRDMs, avgCueRDM] = compute_rdms(data.cueId, @(x1, x2) x1 ~= x2, data, metadata, which_trials);
[contextRDMs, avgContextRDM] = compute_rdms(data.contextId, @(x1, x2) x1 ~= x2, data, metadata, which_trials);
[bothRDMs, avgBothRDM] = compute_rdms(data.contextId * 10 + data.cueId, @(x1, x2) x1 ~= x2, data, metadata, which_trials);

ttest_means = [];
ttest_sems = [];
ttest_ps = [];
ttest_ts = [];

all_fisher_rhos = {};

for neural_idx = 1:numel(Neural) % for each ROI
    fprintf('\n---------- at ROI %s -------------------\n', Neural(neural_idx).name);

    fisher_rhos = []; % for each subject, does neural activity predict test choices, under the different structures?

    for subj = 1:metadata.N % for each subject
        subject = metadata.subjects(subj); % 'con001' ... 'con025'
       
        subj_rhos = []; % M1,M2,M1' rho's (similarity between neural and M1/M2/M1' RDMs) for all runs, i.e. 27 values
        subj_odds = []; % M1,M2,M1' posterior odds based on test choices for all runs, i.e. 27 values

        for run = 1:metadata.runsPerSubject % for each run
            fprintf('subject %d, run %d\n', subj, run);

            t1_mask = data.runId(t1) == run;
            t2_mask = data.runId(t2) == run;
            run_mask = t1_mask & t2_mask & t1 > t2 & t2 > 5; % skip first 5 trials
        
            % how similar is the run RDM (in this ROI) to a M1/M2/M4 RDM (based on attention to cue, context, or both)
            %
            neural_RDM = Neural(neural_idx).RDMs(:,:,subj);

            M1_RDM = cueRDMs(:,:,subj);
            M2_RDM = bothRDMs(:,:,subj);
            M4_RDM = contextRDMs(:,:,subj); % technically M1' 

            neural_subRDM = neural_RDM(run_mask);
            M1_subRDM = M1_RDM(run_mask);
            M2_subRDM = M2_RDM(run_mask);
            M4_subRDM = M4_RDM(run_mask);
            
            rho_M1 = corr(neural_subRDM, M1_subRDM, 'type', 'Spearman');
            rho_M2 = corr(neural_subRDM, M2_subRDM, 'type', 'Spearman');
            rho_M4 = corr(neural_subRDM, M4_subRDM, 'type', 'Spearman');
            assert(~isnan(rho_M1));
            assert(~isnan(rho_M2));
            assert(~isnan(rho_M4));

            M1_rhos(subj, run, neural_idx, 1) = rho_M1; 
            M2_rhos(subj, run, neural_idx, 1) = rho_M2; 
            M4_rhos(subj, run, neural_idx, 1) = rho_M4; 
            subj_rhos = [subj_rhos, rho_M1, rho_M2, rho_M4];

            % simulate run under each structure alone
            %
            which_run = data.which_rows & strcmp(data.participant, subject) & data.runId == run;
            assert(sum(which_run) == metadata.trialsPerRun);

            M1_simulated = simulate_subjects(data, metadata, params, [1 0 0 0 0], which_run, false); % note we use params optimized for all structures together
            M2_simulated = simulate_subjects(data, metadata, params, [0 1 0 0 0], which_run, false);
            M4_simulated = simulate_subjects(data, metadata, params, [0 0 0 1 0], which_run, false);

            % get test choice likelihood under each structure 
            %
            which_run_test = data.which_rows & ~data.isTrain & ~data.timeout & strcmp(data.participant, subject) & data.runId == run;
            assert(sum(which_run_test) <= metadata.testTrialsPerRun);

            sims = {M1_simulated, M2_simulated, M4_simulated};

            struct_logliks = []; % test choice log likelihood under each structure, i.e. log P(test choices | M)
            for sim = sims % for each structure
                sim = sim{1};

                p = sim.pred(which_run_test);
                X = data.chose_sick(which_run_test);
                test_liks = binopdf(X, 1, p);
                assert(numel(test_liks) == numel(X));

                avg_loglik = mean(log(test_liks));
                struct_logliks = [struct_logliks, avg_loglik];
            end

            struct_posts = exp(struct_logliks) / sum(exp(struct_logliks)); % posterior over structures, i.e. P(M | test choices), assuming uniform prior P(M)
            struct_odds = log(struct_posts ./ (1 - struct_posts)); % convert to log odds to make [-inf, +inf] for correlation 

            subj_odds = [subj_odds, struct_odds];
        end

        assert(numel(subj_rhos) == 27);
        assert(numel(subj_odds) == 27);

        % correlate rho's and odds 
        % => does the similarity between structure RDMs and neural RDMs correlate with 
        % the probability of the structures, according to the subject's test choices?
        %
        rho = corr(subj_rhos', subj_odds', 'type', 'Spearman'); % TODO Pearson?
        fisher_rho = atanh(rho);
        fisher_rhos = [fisher_rhos, fisher_rho]; % collect for all subjects

        fprintf('    rho = %f\n', rho);
    end

    % t-test the rho's against 0 for all subjects
    % => does neural activity predict test choices, under the different structures, across subjects?
    %
    [h, p, ci, stats] = ttest(fisher_rhos);
    ttest_ts = [ttest_ts; stats.tstat];
    ttest_ps = [ttest_ps; p];
    ttest_means = [ttest_means; mean(fisher_rhos)];
    ttest_sems = [ttest_sems; (ci(2) - ci(1)) / 2];
    all_fisher_rhos{neural_idx} = fisher_rhos;

    fprintf('-> stats for ROI %s:\n', Neural(neural_idx).name);
    fisher_rhos
    p 
    stats
end

