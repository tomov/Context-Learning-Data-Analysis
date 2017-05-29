function [results, results_options] = fit_params()

% Fit the hyperparameters of the model based on behavioral data
%
% OUTPUT: 
% results = struct array with resulting free parameters for each version of
%           the model that was fit
% results_options = options for the particular version of the model, e.g.
%           was it random effects vs. fixed effects, or was it fmri data or
%           pilot data

results = struct([]); % resulting parameter fits
results_options = struct([]); % settings for each version of the model, e.g. which data we used this time
results_ord = 0; % ordinal = which version of the model we're fitting now

% fit behavioral (0) or fMRI (1) data
%
for isFmriData = [0 1]
    
    options = struct; % options for current result
    options.isFmriData = isFmriData;
    
    if isFmriData
        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
    else
        [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
    end    
    
    % fixed effects (1) or random effects (0)
    % fixed effects = have one "supersubject" and fit one set of parameters for her, 
    % random effects = fit each individual subject separately with her own set of parameters
    %
    for fixedEffects = [1 0]
        
        options.fixedEffects = fixedEffects;
        
        % what we pass to mfit_optimize telling it which rows to look at
        % when fitting each set of parameters. In our case, each set of
        % rows corresponds to the trials for a particular subject, or for
        % all subjects
        %
        mfit_data = struct([]);
        
        if fixedEffects
            % for fixed effects, pass a binary mask that says "use all
            % rows". This means that we'll be fitting each set of
            % parameters to all subjects at the same time
            %
            mfit_data(1).which_rows = data.which_rows;
            % # of data points
            mfit_data(1).N = sum(mfit_data(1).which_rows);
        else
            % for random effects, pass N binary masks (N = # subjects),
            % with each one corresponding to the rows for that subject.
            % This means we'll be fitting a separate set of parameters for
            % each subject
            %
            subj_idx = 0;
            for who = metadata.subjects
                subj_idx = subj_idx + 1;
                mfit_data(subj_idx).which_rows = strcmp(data.participant, who);
                % # of data points
                mfit_data(subj_idx).N = sum(mfit_data(subj_idx).which_rows);
            end
        end
        
        % hypothesis spaces of causal sturctures to fit
        %
        which_structuress = {[1 1 1 0], [1 0 0 0], [0 1 0 0], [0 0 1 0], [1 1 0 0], [1 0 1 0], [0 1 1 0]}; 
        
        for which_structures = which_structuress
            which_structures = which_structures{1};
            
            options.which_structures = which_structures;
          
            %
            % Fit model
            %

            % create parameter structure
            %
            param = [];

            param(1).name = 'prior variance';
            param(1).logpdf = @(x) 1;  % log density function for prior
            param(1).lb = 0; % lower bound
            param(1).ub = 1; % upper bound TODO more?

            param(2).name = 'inverse softmax temperature'; 
            param(2).logpdf = @(x) 1;  % log density function for prior
            param(2).lb = 0;
            param(2).ub = 10; % can't make it too large b/c you get prob = 0 and -Inf likelihiood which fmincon doesn't like

            % P(data | model)
            % the likelihood function takes in all the subject data every
            % time (data and metadata); mfit_data tells it which subjects
            % to simulate by passing a bitmask which includes only rows for
            % this subject / those subjects
            %
            likfun = @(params, mfit_data) model_likfun(data, metadata, params, which_structures, mfit_data.which_rows);
            
            % run optimization
            %
            nstarts = 5;    % number of random parameter initializations 
            results_ord = results_ord + 1;
            fprintf('\nFitting #%d with options:\n', results_ord);
            disp(options);
            
            results(results_ord) = mfit_optimize(likfun, param, mfit_data, nstarts);
            results_options(results_ord) = options;

            fprintf('  BIC = %.4f, AIC = %.4f, params = \n', results(results_ord).bic, results(results_ord).aic);
            disp(results(results_ord).x);
        end
    end
end
