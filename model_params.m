function [params, which_structures] = model_params(params_file, params_idx)
% Load the parameters that we use to simulate subjects.
% Result passed to simulate_subjects()

if nargin < 2 || isempty(params_idx)
    params_idx = 1;
end

load(params_file, 'results', 'results_options');
fprintf('using params from file %s\n', params_file);
params = results(params_idx).x;
options = results_options(params_idx);
which_structures = options.which_structures;

disp('Which structures:');
disp(which_structures);
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);

