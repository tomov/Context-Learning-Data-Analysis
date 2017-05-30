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
        
        %plot_behavior(data, metadata);
        %plot_behavior(data, metadata, [0.1249 2.0064], [1 1 1 0]);
        
        % Plot for a sepecific set of parameters -- the random effects
        % parameters fit to the fMRI subjects behavior
        %
        load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');
        result = results(22);
        options = results_options(22);
        assert(options.isFmriData == true);
        assert(~options.fixedEffects);
        assert(isequal(options.which_structures, [1 1 1 0]));
        
        plot_behavior(data, metadata, result.x, options.which_structures);

        
    case 'parameter_fits'
        % Compare different hyperparameters
        %
        
        fit_params_summary_filename = fullfile('results', 'fit_params_summary.csv');
        if exist(fit_params_summary_filename, 'file') ~= 2
            T = plot_parameterFits();
            writetable(T, fit_params_summary_filename, 'WriteRowNames', true);
        else
            fprintf('Loading summary from file %s\n', fit_params_summary_filename);
            T = readtable(fit_params_summary_filename, 'ReadRowNames', true);
        end        
        disp(T);

    otherwise
        assert(false, 'invalid plotname -- should be one of the above');
        
end


