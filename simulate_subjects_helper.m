function [data, metadata, simulated, params, options, results, results_options] = simulate_subjects_helper(isFmri, params_file, params_idx, which_structures, which_rows)

% helper for running model simulations
%
if nargin < 1
    isFmri = true;
end
if nargin < 2 || isempty(params_file)
    %params_file = fullfile('results', 'fit_params_results.mat');
    params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
end
if nargin < 3 || isempty(params_idx)
    params_idx = 1;
end
if nargin < 4 || isempty(which_structures)
    %which_structures = [1 1 1 0]; % M1, M2, M3
    which_structures = logical([1 1 0 1 0]); % M1, M2, M1'
end

% Load data
%
if isFmri
    % fmri subjects
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
    disp('=== simulate subjects helper, fMRI csv\n');
else
    % pilot subjects
    [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
    disp('=== simulate subjects helper, pilot csv\n');
end

if nargin < 5 || isempty(which_rows)
    which_rows = data.which_rows;
end

% Load parameters
%
%load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_5_prior.mat'), 'results', 'results_options');
load(params_file, 'results', 'results_options');
fprintf('using params from file %s\n', params_file);
params = results(params_idx).x;
options = results_options(params_idx);

%params = [0.1249 2.0064];

disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);

% safeguards
%assert(options.isFmriData == true);
%assert(~options.fixedEffects);
%assert(isequal(options.which_structures, which_structures)); % which_structures provided for sanity check only

% Run the model with the parameters
%
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows);

