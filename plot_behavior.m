function plot_behavior(data, metadata, params, which_structures)

% Plot the subject test trial behavior, learning curve, model posteriors,
% and a few other goodies. These are the main behavioral plots.
% 
% INPUT 
% data, metadata = subject data and metadata as output by load_data
% params = vector of hyperparameters:
%    params(:, 1) = prior_variance
%    params(:, 2) = inv_softmax_temp
%          for fixed effects, only have 1 row and these parameters will be
%          used for all subjects. For random effects, have a separate row
%          with the parameters for each subject
% which_structures = which causal structures to use, as a logical vector,
%                    e.g. [1 1 1 0] = M1, M2, M3
%

if nargin < 4 || isempty(which_structures)
    which_structures = [1 1 1 0]; % by defulat, use M1 M2 M3
end

if nargin < 3 || isempty(params)
    params = [0.1249 2.0064]; % by default, use the params from the pilot fit
end



utils; % include some nifty lambdas

% First simulate the subjects with the causal structure model
%
simulated = simulate_subjects(data, metadata, params, which_structures);

next_subplot_idx = 1; % so you can reorder them by simply rearranging the code

%
% Outcome probabilities in training phase
%

Ms = [];
SEMs = [];

for context = metadata.contextRoles
    which = data.which_rows & data.isTrain == 1 & strcmp(data.contextRole, context);

    x1c1 = strcmp(data.corrAns(which & data.cueId == 0 & data.contextId == 0), 'left');
    x1c2 = strcmp(data.corrAns(which & data.cueId == 0 & data.contextId == 1), 'left');
    x2c1 = strcmp(data.corrAns(which & data.cueId == 1 & data.contextId == 0), 'left');
    x2c2 = strcmp(data.corrAns(which & data.cueId == 1 & data.contextId == 1), 'left');

    M = get_means(x1c1, x1c2, x2c1, x2c2);
    SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
    %M = mean([x1c1 x1c2 x2c1 x2c2]);
    %SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
    Ms = [Ms; M];
    SEMs = [SEMs; SEM];
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(Ms, SEMs, 1, metadata.contextRoles, 'P(sick outcome) in training');
ylabel('Sick probability');
legend({'x_1c_1', 'x_1c_2', 'x_2c_1', 'x_2c_2'});


%
% Choice probabilities in test phase for SUBJECTS
% This is the final figure we care about
%

Ms = [];
SEMs = [];
for context = metadata.contextRoles
    which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

    x1c1 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 0), 'left');
    x1c2 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 2), 'left');
    x2c1 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 0), 'left');
    x2c2 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 2), 'left');

%    M = mean([x1c1 x1c2 x2c1 x2c2]);
%    SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
    M = get_means(x1c1, x1c2, x2c1, x2c2);
    SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
    Ms = [Ms; M];
    SEMs = [SEMs; SEM];
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(Ms, SEMs, 1, metadata.contextRoles, 'Subject P(choose sick) in test');
ylabel('Sick probability');
legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'});


%
% Choice probabilities in test phase for MODEL (based on the actual
% decisions the model made)
% This is the final figure we care about
%

Ms = [];
SEMs = [];
for context = metadata.contextRoles
    which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);


    x1c1 = strcmp(simulated.keys(which & data.cueId == 0 & data.contextId == 0), 'left');
    x1c2 = strcmp(simulated.keys(which & data.cueId == 0 & data.contextId == 2), 'left');
    x2c1 = strcmp(simulated.keys(which & data.cueId == 2 & data.contextId == 0), 'left');
    x2c2 = strcmp(simulated.keys(which & data.cueId == 2 & data.contextId == 2), 'left');

    %M = mean([x1c1 x1c2 x2c1 x2c2]);
    %SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
    M = get_means(x1c1, x1c2, x2c1, x2c2);
    SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
    Ms = [Ms; M];
    SEMs = [SEMs; SEM];
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(Ms, SEMs, 1, metadata.contextRoles, 'Model P(choose sick) in test');
ylabel('Sick probability');
legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'});


%
% TRUE Choice probabilities in test phase for MODEL
%

Ms = [];
SEMs = [];
for context = metadata.contextRoles
    which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

    x1c1 = simulated.pred(which & data.cueId == 0 & data.contextId == 0);
    x1c2 = simulated.pred(which & data.cueId == 0 & data.contextId == 2);
    x2c1 = simulated.pred(which & data.cueId == 2 & data.contextId == 0);
    x2c2 = simulated.pred(which & data.cueId == 2 & data.contextId == 2);

    %M = mean([x1c1 x1c2 x2c1 x2c2]);
    %SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
    M = get_means(x1c1, x1c2, x2c1, x2c2);
    SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
    Ms = [Ms; M];
    SEMs = [SEMs; SEM];
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(Ms, SEMs, 1, metadata.contextRoles, 'Ideal model P(choose sick) in test');
ylabel('Sick probability');
legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'});




%
% Per-trial accuracy across all subjects & runs (training only)
% compared against the model
%

human_correct = [];
model_correct = [];

for n = 1:metadata.trainingTrialsPerRun
    which = data.which_rows & data.isTrain & data.trialId == n;

    human_corr_n = strcmp(data.response.keys(which), data.corrAns(which));
    model_corr_n = strcmp(simulated.keys(which), data.corrAns(which));
    human_correct = [human_correct mean(human_corr_n)];
    model_correct = [model_correct mean(model_corr_n)];
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
plot(model_correct, 'o-', 'LineWidth', 2); % == mean(human_correct_all_runs)
hold on;
plot(human_correct, 'o-', 'LineWidth', 2); % == mean(model_correct_all_runs)
hold off;
legend({'model', 'subject'});
title('Per-trial accuracy');
xlabel('trial #');
ylabel('accuracy');


%
% Model per-trial posterior probability P(M | ...) for each condition
%

for condition = metadata.contextRoles

    P = [];
    for n = 1:metadata.trainingTrialsPerRun
        which = data.which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

        P1_n = simulated.P1(which);
        P2_n = simulated.P2(which);
        P3_n = simulated.P3(which);
        P4_n = simulated.P4(which);

        P_n = [];  % only include the posteriors from models we care about
        if which_structures(1), P_n = [P_n mean(P1_n)]; end
        if which_structures(2), P_n = [P_n mean(P2_n)]; end
        if which_structures(3), P_n = [P_n mean(P3_n)]; end
        if which_structures(4), P_n = [P_n mean(P4_n)]; end
        P = [P; P_n];
    end

    subplot(3, 5, next_subplot_idx);
    next_subplot_idx = next_subplot_idx + 1;
    plot(P, 'o-', 'LineWidth', 2);
    xlabel('n (trial #)');
    ylabel('P(M | h_{1:n})');
    title(strcat('Posterior after each trial for ', {' '}, condition));

    Ms = {}; % only include models we care about
    if which_structures(1), Ms = [Ms, {'M1'}]; end
    if which_structures(2), Ms = [Ms, {'M2'}]; end
    if which_structures(3), Ms = [Ms, {'M3'}]; end
    if which_structures(4), Ms = [Ms, {'M4'}]; end
    legend(Ms);

end


%
% Model per-trial weight matrix ww for each condition
%
for condition = metadata.contextRoles

    ww = [];
    for n = 1:metadata.trainingTrialsPerRun
        which = data.which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

        ww1_n = simulated.ww1(which, :);
        ww2_n = simulated.ww2(which, :);
        ww3_n = simulated.ww3(which, :);
        ww4_n = simulated.ww4(which, :);

        ww_n = []; % only include the weights from models we care about
        if which_structures(1), ww_n = [ww_n mean(ww1_n)]; end
        if which_structures(2), ww_n = [ww_n mean(ww2_n)]; end
        if which_structures(3), ww_n = [ww_n mean(ww3_n)]; end
        if which_structures(4), ww_n = [ww_n mean(ww4_n)]; end
        ww = [ww; ww_n];
    end

    subplot(3, 5, next_subplot_idx);
    next_subplot_idx = next_subplot_idx + 1;
    plot(ww, 'o-', 'LineWidth', 1);
    xlabel('n (trial #)');
    ylabel('ww_n');
    title(strcat('Weights after each trial for ', {' '}, condition));

    Ms = {}; % only include models we care about
    if which_structures(1), Ms = [Ms, {'M1, x1'}, {'M1, x2'}]; end
    if which_structures(2), Ms = [Ms, {'M2, x1c1'}, {'M2, x2c1'}, {'M2, x1c2'}, {'M2, x2c2'}]; end
    if which_structures(3), Ms = [Ms, {'M3, x1'}, {'M3, x2'}, {'M3, c1'}, {'M3, c2'}]; end
    if which_structures(4), Ms = [Ms, {'M4, c1'}, {'M4, c2'}]; end
    legend(Ms);

end

%
% D_KL == bayesian surprise on each training trial, by whether trial was
% correct
%
P = [simulated.P1 simulated.P2 simulated.P3 simulated.P4];
Q = [simulated.Q1 simulated.Q2 simulated.Q3 simulated.Q4];

means = zeros(1,2);
sems = zeros(1,2);
sem = @(x) std(x) / sqrt(length(x));
for correct = 1:-1:0
    which = data.which_rows & data.isTrain & data.response.corr == correct;

    means(1, 2 - correct) = mean(simulated.surprise(which)); 
    sems(1, 2 - correct) = sem(simulated.surprise(which));
end

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(means, sems);
legend({'correct', 'wrong'});
ylabel('D_{KL} on current trial');

%
% D_KL == bayesian surprise on each training trial, by whether next trial
% was correct
%

which_next_trials = data.which_rows & data.isTrain & data.trialId ~= 1;
which_prev_trials = data.which_rows & data.isTrain & data.trialId ~= 20;
prev_trials_surprise = simulated.surprise(which_prev_trials);
next_trials_corr = data.response.corr(which_next_trials);


means = [mean(prev_trials_surprise(next_trials_corr == 1)) mean(prev_trials_surprise(next_trials_corr == 0))];
sems = [sem(prev_trials_surprise(next_trials_corr == 1)) sem(prev_trials_surprise(next_trials_corr == 0))];

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(means, sems);
legend({'correct', 'wrong'});
ylabel('D_{KL} on previous trial');




%
% value in sick vs. not sick
%
%{
outcomes = strcmp(sick(data.which_rows & data.isTrain), 'Yes');
values = simulated.values(data.which_rows & data.isTrain);

means = [mean(values(outcomes == 1)) mean(values(outcomes == 0))];
sems = [sem(values(outcomes == 1)) sem(values(outcomes == 0))];

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(means, sems);
legend({'sick', 'not sick'});
ylabel('value');
%}


%
% value in chose sick vs. chose not sick, on wrong trials
%

wrong_trials = data.which_rows & data.isTrain & data.response.corr == 0;
subject_choices = strcmp(data.response.keys(wrong_trials), 'left');
values = simulated.values(wrong_trials);

means = [mean(values(subject_choices == 1)) mean(values(subject_choices == 0))];
sems = [sem(values(subject_choices == 1)) sem(values(subject_choices == 0))];

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(means, sems);
title('Wrong Trials');
legend({'chose sick', 'chose not sick'});
ylabel('value');

%
% value in chose sick vs. chose not sick, on correct trials
%

correct_trials = data.which_rows & data.isTrain & data.response.corr == 1;
subject_choices = strcmp(data.response.keys(correct_trials), 'left');
values = simulated.values(correct_trials);

means = [mean(values(subject_choices == 1)) mean(values(subject_choices == 0))];
sems = [sem(values(subject_choices == 1)) sem(values(subject_choices == 0))];

subplot(3, 5, next_subplot_idx);
next_subplot_idx = next_subplot_idx + 1;
barweb(means, sems);
title('Correct Trials');
legend({'chose sick', 'chose not sick'});
ylabel('value');
