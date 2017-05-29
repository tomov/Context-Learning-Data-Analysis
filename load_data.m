function [data, metadata] = load_data(filename, isFmriData, goodSubjects_ords)

% Load the behavioral data from the .csv file generated by
% snippets/parse.py from the psychopy code.
% It's assumed to be in "wide format" csv for all subjects.
%
% Example USAGE: 
% [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects())
% [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false)
%
% INPUT:
% filename = path to csv file e.g. 'data/fmri.csv' or 'data/pilot.csv'
% isFmriData = boolean flag whether this is data from the fmri sessions
%                or the behavioral pilot sessions. Those have different
%                formats so it's important to get this right.
% goodSubject_ords (optional) = vector of indices of the subjects to include in the
%                analysis. If empty, assumes all subjects are included.
%                Indices correspond to the ordinal of the subject in the
%                csv file
%
% OUTPUT:
% metadata = struct with data about the data
% data = struct with fields where each field is a column vector; each row
%        of the column vectors corresponds to a single trial of a given
%        subject. Then you select particular trials using logical arrays
%        based. Some of the more important ones are:
%    .which_rows = which rows to include in the analysis. By default all of
%            them. Every logical vector must be &'d with which_rows
%    .chose_sick = whether the subject chose 'Sick' on that trial
%    .timeout = whether the subject timed out
%    .isTrain = is this a training trial (vs. a test trial)
%    .runId = which run/round/block of the session
%    .trialId = trial ordinal; test trials start over (1..4)
%    .contextRole = condition
%    .participant = subject id as entered in psychopy
%    .corrAns = the correct answer on that trial
%    .sick = 'Yes' or 'No', whether the outcome was sick
%    .response.corr = was the subject correct (training trials only)
%    .response.keys = 'left' if subject said 'sick', 'right' if said 'not
%                     sick', 'None' if subject timed out
%

% Set default arguments
%
if nargin < 3
    goodSubjects_ords = [];
end

% Load from CSV file
%
if isFmriData
    % behavioral data from fMRI subjects
    %
    format = '%d %s %s %s %d %s %s %d %d %s %s %s %f %d %s %s %d %d %d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
    [data.drop, data.participant, data.session, data.mriMode, data.isPractice, data.runFilename,...
        data.contextRole, data.contextId, data.cueId, data.sick, data.corrAns, ...
        data.response.keys, data.response.rt, data.response.corr, data.restaurant, data.food, ...
        data.isTrain, data.runId, data.trialId, ...
        data.trialStartWallTime, data.actualChoiceOnset, data.actualChoiceOffset, ...
        data.actualIsiOnset, data.actualIsiOffset, data.actualFeedbackOnset, data.actualFeedbackOffset, ...
        data.actualItiOnset, data.actualItiOffset, data.actualItiDuration, data.itiDriftAdjustment, ...
        data.trialEndWallTime, data.stimOnset, data.stimOffset, data.itiDuration, data.itiOffset] = ...
        textread(filename, format, 'delimiter', ',', 'headerlines', 1);
else
    % behavioral data from pilot subjects
    %
    format = '%s %s %s %d %s %s %s %d %d %s %s %s %f %d %s %s %d %d %d';

    [data.participant, data.session, data.mriMode, data.isPractice, data.restaurantsReshuffled, data.foodsReshuffled, data.contextRole, data.contextId, data.cueId, data.sick, data.corrAns, data.response.keys, data.response.rt, data.response.corr, data.restaurant, data.food, data.isTrain, data.runId, data.trialId] = ...
        textread(filename, format, 'delimiter', ',', 'headerlines', 1);
end
    
data.which_rows = logical(true(size(data.participant))); % by default, include all rows
data.chose_sick = strcmp(data.response.keys, 'left');
data.timeout = strcmp(data.response.keys, 'None');
data.isTest = ~data.isTrain;
data.outcome = strcmp(data.sick, 'Yes');
data.condition = data.contextRole;

metadata.allSubjects = unique(data.participant)';

% Only include the "good" subjects
%
if ~isempty(goodSubjects_ords)
    goodSubjects = metadata.allSubjects(goodSubjects_ords);
    data.which_rows = data.which_rows & ismember(data.participant, goodSubjects);
end

if isFmriData
    data.no_response = cellfun(@isempty, data.actualChoiceOffset);
    assert(isequal(data.timeout, data.no_response), 'these should be one and the same -- maybe bug in psychopy code');

    % we shouldn't be including any dropped trials in the analysis
    assert(sum(data.drop(data.which_rows)) == 0, 'bad trials found -- perhaps need to pass a list of good subjects');
else
    data.no_response = data.timeout;
end

% Set the metadata
%
metadata.runsPerContext = 3; % = blocks per context = runs per context = runs / 3
metadata.runsPerSubject = 9;
metadata.trialRepsPerRun = 5; % = training trials per run / 4
metadata.trainingTrialsPerRun = 20;
metadata.testTrialsPerRun = 4;
metadata.trialsPerRun = metadata.trainingTrialsPerRun + metadata.testTrialsPerRun;
metadata.subjects = unique(data.participant(data.which_rows))'; % only the included subjects
metadata.N = numel(metadata.subjects);
metadata.contextRoles = {'irrelevant', 'modulatory', 'additive'}; % in order
metadata.conditions = metadata.contextRoles;
assert(isequal(sort(metadata.contextRoles), sort(unique(data.contextRole))'), 'maybe extra contextRoles in data');

% Set some more data
%
data.newTrialId = data.trialId + data.isTest * metadata.trainingTrialsPerRun; % numbered 1..24 for the run (otherwise it's 1..20 for training trials and 1..4 for test trials)

