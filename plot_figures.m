function plot_figures(plotname, isFmriData)

% set default params
%
if nargin < 2 || isempty(isFmriData)
    isFmriData = true;
end

% load behavioral data
%
if isFmriData
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
else
    [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
end

% Plot figure according to plotname
%
switch plotname
    
    case 'sanity'
        % Basic plot to make sure the design was correct and subjects
        % understood the task
        %
        
        figure;
        plot_sanity(data, metadata);
        
    case 'behavior'
        % The main behavioral plots -- test choices, model choices,
        % accuracies, etc
        %
        
        figure;
        plot_behavior(data, metadata);
        %plot_behavior(data, metadata, [0.1249 2.0064], [1 1 1 0]);
        
    case 'parameter_fits'
        % Compare different hyperparameters
        %
   
        % Check if we already computed the parameter fits.
        % If we haven't, compute them (using mfit) and save them.
        %
        fit_params_filename = fullfile('results', 'fit_params_results.mat');
        if exist(fit_params_filename, 'file') ~= 2
            fprintf('Could not find saved fit param results in %s; recomputing...\n', fit_params_filename);
            [results, results_options] = fit_params();
            save(fit_params_filename, 'results', 'results_options');
        else
            fprintf('Loading fit param results from %s...\n', fit_params_filename);
            load(fit_params_filename, 'results', 'results_options');
        end

        
        table_rowNames = cell(numel(results), 1);
        table_colNames = {'BIC', 'loglik', 'num_params', 'corr', 'corr_test', 'row'};
        table = nan(numel(results), numel(table_colNames));

        struct_names = {'M1, ', 'M2, ', 'M3, ', 'M4, '};
        pilot_or_fmri = {'pilot', 'fmri'};
        fixed_or_random_effects = {'random effects', 'fixed effects'};

        assert(numel(results) == numel(results_options));
        for i = 1:numel(results)
            isFmriData = results_options(i).isFmriData;
            fixedEffects = results_options(i).fixedEffects;
            which_structures = results_options(i).which_structures;
            
            rowName = sprintf('%s; %s; %s', pilot_or_fmri{isFmriData + 1}, ...
                fixed_or_random_effects{fixedEffects + 1}, strcat(struct_names{logical(which_structures)}));
            table_rowNames(i) = {rowName};
            
            result = results(i);
            table(i, 1) = sum(result.bic);
            table(i, 2) = sum(result.loglik);
            table(i, 3) = numel(result.x);
            table(i, end) = i;
        end

        T = array2table(table, 'RowNames', table_rowNames, 'VariableNames', table_colNames);
        disp(T);

    otherwise
        assert(false, 'invalid plotname -- should be one of the above');
        
end


