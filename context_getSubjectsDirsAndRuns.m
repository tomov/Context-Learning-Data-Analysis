function [ subjects, subjdirs, nRuns ] = context_getSubjectsDirsAndRuns()

% Get the list of subjects, subject directories, and number of runs for the
% fMRI GLM code
%

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% the participant id as entered in psychopy
subjects = {'con001', 'con002', 'con003', 'con004', ...
            'con005', 'con006', 'con007', 'con008', ...
            'con009', 'con010', 'con011', 'con012', ...
            'con013', 'con014', 'con015', 'con017', ...
            'con018', 'con021', 'con022', 'con023', ...
            'con024', 'con025', 'con026', 'con027', ...
            'con028'}; 
        
% should be identical to the list of subjects in the csv file
% and in the same order
% this is a basic assumption for getGoodSubjects() to work
% we are listing them here explicitly as a sanity check for the csv file
%
assert(mean(strcmp(subjects, unique(data.participant)')) == 1);

% the names of the CORRESPONDING directories from CBS central
subjdirs = {'161030_con001', '161030_con002', '161101_CON_003', '161106_CON_004', ...
            '161106_CON_005', '161106_CON_006', '161106_CON_007', '161108_CON_008', ...
            '161112_CON_009', '161112_CON_010', '161112_CON_011', '161112_CON_012', ...
            '161203_CON_013', '161203_CON_014', '161203_CON_015', '161203_CON_017', ...
            '161204_CON_018', '161204_CON_021', '161204_CON_022', '161204_CON_023', ...
            '161204_CON_024', '161205_CON_025', '161206_CON_026', '161210_CON_027', ...
	    '161211_CON_028'};

% assumes runs are always in order: 1,2,3,4,...
nRuns = {9, 9, 5, 9, ...
         9, 9, 9, 9, ...
         6, 9, 9, 9, ...
         9, 9, 9, 9, ...
         9, 9, 9, 9, ...
         9, 9, 9, 9, ...
	 9}; % runs per subject



