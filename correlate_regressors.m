% Correlate different potential regressors
%

% Load data & simulate subjects
%
%{
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

load(fullfile('results', 'fit_params_results_BAD.mat'), 'results', 'results_options');
result = results(22);
options = results_options(22);
assert(options.isFmriData == true);
assert(~options.fixedEffects);
assert(isequal(options.which_structures, [1 1 1 0]));

simulated = simulate_subjects(data, metadata, result.x, options.which_structures);        
%}

%
% Correlate lambdas i.e. 1/precision for each causal structure
%

which_rows = data.which_rows & data.isTrain;
lambdas = simulated.lambdas(:, logical(options.which_structures));

% For all subjects
%
[r_lambdas, p_lambdas] = corrcoef(lambdas(which_rows, :))

% For each subject separately
%
rs_lambdas = nan(size(lambdas, 2), size(lambdas, 2), metadata.N);
s_id = 0;
for who = metadata.subjects
    s_id = s_id + 1;
    r = corrcoef(lambdas(which_rows & strcmp(data.participant, who), :));
    rs_lambdas(:, :, s_id) = r;
end


%
% Correlate posteriors
%

which_rows = data.which_rows & data.isTrain;
P = simulated.P(:, logical(options.which_structures));

% For all subjects
%
[r_P, p_P] = corrcoef(P(which_rows, :))

% For each subject separately
%
rs_P = nan(size(P, 2), size(P, 2), metadata.N);
s_id = 0;
for who = metadata.subjects
    s_id = s_id + 1;
    r = corrcoef(P(which_rows & strcmp(data.participant, who), :));
    rs_P(:, :, s_id) = r;
end

