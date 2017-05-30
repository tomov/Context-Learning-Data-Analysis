function T = plot_parameterFits()

% Show hyperparameter fits for different sets of data, differnet settings
% (fixed vs. random effects) and different hypothesis spaces.
%
% OUTPUT:
% A table where each row represents a different set of fit parameters and
% each column is some attribute characterizing it (BIC, etc)
%

% Check if we already computed the parameter fits.
% If we haven't, compute them (using mfit) and save them.
%
fit_params_filename = fullfile('results', 'fit_params_results.mat');
if exist(fit_params_filename, 'file') ~= 2
    fprintf('Could not find saved fit param results in %s; recomputing...\n', fit_params_filename);
    [results, results_options, mfit_datas] = fit_params();
    save(fit_params_filename, 'results', 'results_options', 'mfit_datas', '-v7.3');
else
    fprintf('Loading fit param results from %s...\n', fit_params_filename);
    load(fit_params_filename, 'results', 'results_options', 'mfit_datas');
end

% Plot them in a table
%
table_rowNames = cell(numel(results), 1);
table_colNames = {'BIC', 'loglik', 'loglik_pilot_all', 'loglik_pilot_test',  'r_pilot_all', 'p_pilot_all', 'r_pilot_test', 'p_pilot_test',  'loglik_fmri_all', 'loglik_fmri_test', 'r_fmri_all', 'p_fmri_all', 'r_fmri_test', 'p_fmri_test', 'num_params', 'row'};
table = nan(numel(results), numel(table_colNames));

struct_names = {'M1, ', 'M2, ', 'M3, ', 'M4, '};
pilot_or_fmri = {'pilot', 'fmri'};
fixed_or_random_effects = {'random effects', 'fixed effects'};

assert(numel(results) == numel(results_options));
for i = 1:numel(results)
    disp(i);
    
    isFmriData = results_options(i).isFmriData;
    fixedEffects = results_options(i).fixedEffects;
    which_structures = results_options(i).which_structures;

    % Come up with a user-friendly name for each row
    %
    rowName = sprintf('%s; %s; %s', pilot_or_fmri{isFmriData + 1}, ...
        fixed_or_random_effects{fixedEffects + 1}, strcat(struct_names{logical(which_structures)}));
    table_rowNames(i) = {rowName};
    
    % Get the parameter fit attributes
    %
    result = results(i);
    
    % Compute log likelihoods and correlations for pilot/fmri behavioral
    % data vs. model predictions, for all trials and test trials separately
    %
    if isFmriData || fixedEffects
        % Compute for fMRI data. Note that we can't if these are random
        % effects parameters for the pilot data
        %
        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
        simulated = simulate_subjects(data, metadata, result.x, which_structures);
        
        which = data.which_rows & ~data.timeout;
        
        loglik_fmri_all = sum(log(binopdf(data.chose_sick(which), 1, simulated.pred(which))));
        
        loglik_fmri_test = sum(log(binopdf(data.chose_sick(which & ~data.isTrain), 1, simulated.pred(which & ~data.isTrain))));
        
        [r, p] = corrcoef(data.chose_sick(which), simulated.pred(which));
        r_fmri_all = r(1,2);
        p_fmri_all = p(1,2);
        
        [r, p] = corrcoef(data.chose_sick(which & ~data.isTrain), simulated.pred(which & ~data.isTrain));
        r_fmri_test = r(1,2);
        p_fmri_test = p(1,2);
    else
        loglik_fmri_all = NaN;
        loglik_fmri_test = NaN;
        r_fmri_all = NaN;
        p_fmri_test = NaN;
    end
    
    if ~isFmriData || fixedEffects
        % Compute for pilot data. Note that we can't if these are
        % parameters for random effects for fmri data
        %
        [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
        simulated = simulate_subjects(data, metadata, result.x, which_structures);
        
        which = data.which_rows & ~data.timeout;
        
        loglik_pilot_all = sum(log(binopdf(data.chose_sick(which), 1, simulated.pred(which))));
        
        loglik_pilot_test = sum(log(binopdf(data.chose_sick(which & ~data.isTrain), 1, simulated.pred(which & ~data.isTrain))));
        
        [r, p] = corrcoef(data.chose_sick(which), simulated.pred(which));
        r_pilot_all = r(1,2);
        p_pilot_all = p(1,2);
        
        [r, p] = corrcoef(data.chose_sick(which & ~data.isTrain), simulated.pred(which & ~data.isTrain));
        r_pilot_test = r(1,2);
        p_pilot_test = p(1,2);
    else
        loglik_pilot_all = NaN;
        loglik_pilot_test = NaN;
        r_pilot_all = NaN;
        p_pilot_test = NaN;
    end
        
    % Populate the summary table
    %
    table(i, 1) = sum(result.bic);
    table(i, 2) = sum(result.loglik);
    table(i, 3) = loglik_pilot_all;
    table(i, 4) = loglik_pilot_test;
    table(i, 5) = r_pilot_all;
    table(i, 6) = p_pilot_all;
    table(i, 7) = r_pilot_test;
    table(i, 8) = p_pilot_test;
    table(i, 9) = loglik_fmri_all;
    table(i, 10) = loglik_fmri_test;
    table(i, 11) = r_fmri_all;
    table(i, 12) = p_fmri_all;
    table(i, 13) = r_fmri_test;
    table(i, 14) = p_fmri_test;
    table(i, end-1) = numel(result.x);
    table(i, end) = i;
end

T = array2table(table, 'RowNames', table_rowNames, 'VariableNames', table_colNames);
