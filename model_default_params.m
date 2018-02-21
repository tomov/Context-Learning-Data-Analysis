function [params, which_structures] = model_default_params()
% Load the default parameters that we use to simulate subjects.
% Result passed to simulate_subjects()
%
% Currently it's M1, M2, M1'
% with with pilot data, fixed effects
%

params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
params_idx = 1;

load(params_file, 'results', 'results_options');
fprintf('using params from file %s\n', params_file);
params = results(params_idx).x;
options = results_options(params_idx);
which_structures = logical(options.which_structures);

disp('Which structures:');
disp(which_structures);
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);

% safeguards
assert(options.isFmriData == false);
assert(options.fixedEffects == true);
assert(isequal(options.which_structures, [1 1 0 1 0]));
