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

% Load the parameters from the mfit_optimize results
%
load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_100_prior.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
%assert(options.isFmriData == true);
%assert(~options.fixedEffects);
assert(isequal(options.which_structures, [1 1 1 0]));
        

% Plot figure according to plotname
%
switch plotname
    
    case 'sanity'
        figure;
        plot_sanity(data, metadata);
        
    case 'behavior'
        figure;        
        plot_behavior(data, metadata, params, options.which_structures);

    case 'posteriors'
        figure;        
        plot_posteriors(data, metadata, params, options.which_structures);
        
    case 'reliabilities'
        figure;
        plot_reliabilities(data, metadata, params, options.which_structures);
        
    case 'surprise'
        figure;
        plot_surprise(data, metadata, params, options.which_structures);
        
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


