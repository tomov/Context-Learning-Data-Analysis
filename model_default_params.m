function [params, which_structures] = model_default_params()
% Load the default parameters that we use to simulate subjects.
% Result passed to simulate_subjects()
%
% Currently it's M1, M2, M1'
% with with pilot data, fixed effects
%

params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
params_idx = 1;

[params, which_structures] = model_params(params_file, params_idx);

% safeguards
assert(options.isFmriData == false);
assert(options.fixedEffects == true);
assert(isequal(options.which_structures, [1 1 0 1 0]));


assert(false, 'DON''T USE'); 

