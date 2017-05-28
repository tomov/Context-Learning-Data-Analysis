function plot_figures(plotname, isFmriData)

% set default params
%
if nargin < 2 || isempty(isFmriData)
    isFmriData = true;
end

% load behavioral data
%
if isFmriData
    [data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
else
    [data, metadata] = load_data('data/pilot.csv', false);
end

% Plot figure according to plotname
%
switch plotname
    
    case 'sanity'
        
        figure;
        plot_sanity(data, metadata);
        
    case 'behavior'
        
        figure;
        plot_behavior(data, metadata);
        %plot_behavior(data, metadata, [0.1249 2.0064], [1 1 1 0]);
        
    case 'parameter_fits'
   
        % Check if we already computed the parameter fits
        %
        if ~exist('fit_params_output.mat', 'file') == 2
            [results, results_options] = fit_params();
            save('fit_params_output.mat');
        else
            load('fit_params_output.mat');
        end
        
        % TODO plot stuff
        
end

