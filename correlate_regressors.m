% Correlate different potential regressors
%

% Load data & simulate subjects (optionally)
%

preloaded_filename = fullfile('temp', 'correlate_regressors.mat');
if exist(preloaded_filename, 'file') ~= 2
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

    load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');
    result = results(22);
    options = results_options(22);
    assert(options.isFmriData == true);
    assert(~options.fixedEffects);
    assert(isequal(options.which_structures, [1 1 1 0]));

    simulated = simulate_subjects(data, metadata, result.x, options.which_structures);        
    
    save(preloaded_filename);
else
    load(preloaded_filename);
end

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

fisher_rs_lambdas = atanh(rs_lambdas);
for i = 1:3
    for j = 1:i-1
        r = squeeze(fisher_rs_lambdas(i, j, :));
        [h, p] = ttest(r);
        fprintf('Within-subject correlation between lambdas for %s and %s: t-test h = %d, p = %f (mean r = %f)\n', ...
            struct_names{i}, struct_names{j}, h, p, mean(rs_lambdas(i, j, :)));
    end
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

struct_names = {'M1', 'M2', 'M3', 'M4'};

fisher_rs_P = atanh(rs_P);
for i = 1:3
    for j = 1:i-1
        r = squeeze(fisher_rs_P(i, j, :));
        [h, p] = ttest(r);
        fprintf('Within-subject correlation between posteriors for %s and %s: t-test h = %d, p = %f (mean r = %f)\n', ...
            struct_names{i}, struct_names{j}, h, p, mean(rs_P(i, j, :)));
    end
end
