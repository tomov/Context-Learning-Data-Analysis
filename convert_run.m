function [train_x, train_k, train_r, test_x, test_k] = convert_run(data, metadata, who, run)

% Convert the stimulus sequence from a run for a given subject into a
% format that's suitable for the model (i.e. that can be passed to
% model_train.m -- see that for more info)
% 
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% who = the subject id as entered in psychopy, e.g. 'con001'
% run = the run number
%
% OUTPUT:
% train_x, train_k, train_r = the first three inputs to model_train.m (cues, contexts &
%                             outcomes)
% test_x, test_k = same but for model_test.m
%

which_trials = data.which_rows & strcmp(data.participant, who) & data.runId == run;
assert(sum(which_trials) == metadata.trialsPerRun);

which_train = which_trials & data.isTrain;
which_test = which_trials & ~data.isTrain;
assert(sum(which_train) == metadata.trainingTrialsPerRun);
assert(sum(which_test) == metadata.testTrialsPerRun);

% training trials
%
cues = data.cueId(which_train);
N = length(cues); % # of trials
D = 3; % # of stimuli
train_x = zeros(N, D);
train_x(sub2ind(size(train_x), 1:N, cues' + 1)) = 1;
train_k = data.contextId(which_train) + 1;
train_r = strcmp(data.sick(which_train), 'Yes');

% test trials
%
test_cues = data.cueId(which_test);
test_N = length(test_cues); % # of trials
D = 3; % # of stimuli
test_x = zeros(test_N, D);
test_x(sub2ind(size(test_x), 1:test_N, test_cues' + 1)) = 1;
test_k = data.contextId(which_test) + 1;


