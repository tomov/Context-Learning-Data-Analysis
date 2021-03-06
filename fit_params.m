function [results, results_options, mfit_datas] = fit_params(isFmriDataRange, fixedEffectsRange, which_structuress, nstarts, outfile, nparams)

% Fit the hyperparameters of the model based on behavioral data
%
% INPUT:
% isFmriDataRange (optional) = fmri or pilot data or both, e.g. [0 1]
% fixedEffectsRange (optional) fixed or random effects or both, e.g. [0 1]
% which_structuress (optional) = cell array of bitmasks for structure
%                                hypothesis spaces to fit, e.g. {[1 1 1 0], [1 0 0 0]}
% nstarts (optional) = how many starts for mfit_optimize for each model
%                      version
% outfile (optional) = where to save the results
% nparams = how many of the parameters to fit (some models have default parameters)
%
% OUTPUT: 
% results = struct array with resulting parameter fits for each version of
%           the model that was fit
% results_options = options for the particular version of the model, e.g.
%           was it random effects vs. fixed effects, or was it fmri data or
%           pilot data
% mfit_datas = the data struct passed to mfit_optimize for each version of
%              the model

% set some default parameters
%
if nargin < 1 || isempty(isFmriDataRange)
    % Do pilot and fMRI subjects
    %
    isFmriDataRange = [0 1];
end
if nargin < 2 || isempty(fixedEffectsRange)
    % Do both fixed and random effects
    %
    fixedEffectsRange = [1 0];
end
if nargin < 3 || isempty(which_structuress)
    % hypothesis spaces of causal sturctures to fit
    %
    which_structuress = {[1 1 1 0], [1 0 0 0], [0 1 0 0], [0 0 1 0], [1 1 0 0], [1 0 1 0], [0 1 1 0]};    
end
if nargin < 4 || isempty(nstarts)
    % number of random parameter initializations 
    %
    nstarts = 5;
end
if nargin < 5 || isempty(outfile)
    outfile = [];
end
if nargin < 6 || isempty(nparams)
    nparams = 0;
end


tic

% fit behavioral (0) or fMRI (1) data
%
results_ord = 0; % ordinal = which version of the model we're fitting now
for isFmriData = isFmriDataRange
    
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
    for fixedEffects = fixedEffectsRange
        
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
            for who = metadata.subjects % only include "good" subjects
                subj_idx = subj_idx + 1;
                mfit_data(subj_idx).which_rows = strcmp(data.participant, who);
                % # of data points
                mfit_data(subj_idx).N = sum(mfit_data(subj_idx).which_rows);
            end
        end
        
        for which_structures = which_structuress
            which_structures = which_structures{1};
            
            options.which_structures = which_structures;
          
            %
            % Fit model
            %

            % create parameter structure
            %
            param = [];

            if isequal(which_structures, 'simple_Q') 
                param(1).name = 'inverse softmax temperature'; 
                param(1).logpdf = @(x) 1;
                param(1).lb = 0;
                param(1).ub = 10;

            elseif isequal(which_structures, 'Q_learning')
                param(1).name = 'learning rate';
                param(1).logpdf = @(x) 1; 
                param(1).lb = 0;
                param(1).ub = 1;

                param(2).name = 'inverse softmax temperature'; 
                param(2).logpdf = @(x) 1;
                param(2).lb = 0;
                param(2).ub = 10; 

                if nparams >= 3
                    param(3).name = 'initial Q-value'; 
                    param(3).logpdf = @(x) 1;
                    param(3).lb = 0;
                    param(3).ub = 1; 
                end

            elseif isequal(which_structures, 'simple_collins') 
                param(1).name = 'learning rate';
                param(1).logpdf = @(x) 1; 
                param(1).lb = 0;
                param(1).ub = 1;

                param(2).name = 'inverse softmax temperature'; 
                param(2).logpdf = @(x) 1;
                param(2).lb = 0;
                param(2).ub = 10; 

                param(3).name = 'concentration parameter'; 
                param(3).logpdf = @(x) 1;
                param(3).lb = 0;
                param(3).ub = 10; 

                if nparams >= 4
                    param(4).name = 'initial Q-value'; 
                    param(4).logpdf = @(x) 1;
                    param(4).lb = 0;
                    param(4).ub = 1; 
                end

            elseif isequal(which_structures, 'flat_collins') 
                param(1).name = 'learning rate';
                param(1).logpdf = @(x) 1; 
                param(1).lb = 0;
                param(1).ub = 1;

                param(2).name = 'inverse softmax temperature'; 
                param(2).logpdf = @(x) 1;
                param(2).lb = 0;
                param(2).ub = 10; 

                if nparams >= 3
                    param(3).name = 'initial Q-value'; 
                    param(3).logpdf = @(x) 1;
                    param(3).lb = 0;
                    param(3).ub = 1; 
                end

            else
                % Sam's model
                assert(~ischar(which_structures));
                param(1).name = 'prior variance';
                param(1).logpdf = @(x) 1;  % log density function for prior
                param(1).lb = 0; % lower bound
                param(1).ub = 10; % upper bound TODO more?

                param(2).name = 'inverse softmax temperature'; 
                param(2).logpdf = @(x) 1;  % log density function for prior
                param(2).lb = 0;
                param(2).ub = 10; % can't make it too large b/c you get prob = 0 and -Inf likelihiood which fmincon doesn't like

                if nparams >= 3
                    param(3).name = 'diffusion_variance'; 
                    param(3).logpdf = @(x) 1;  % log density function for prior
                    param(3).lb = 0;
                    param(3).ub = 1;
                end

                if nparams >= 4
                    param(4).name = 'initial_weight'; 
                    param(4).logpdf = @(x) 1;  % log density function for prior
                    param(4).lb = 0;
                    param(4).ub = 1;
                end
            end


            % P(behavioral data | model)
            % the likelihood function takes in all the subject data every
            % time (data and metadata); mfit_data tells it which subjects
            % to simulate by passing a bitmask which includes only rows for
            % this subject / those subjects
            %
            assert(isequal(mfit_data.which_rows | data.which_rows, data.which_rows)); % make sure we're not including "bad" subjects
            likfun = @(params, mfit_data) model_likfun(data, metadata, params, which_structures, mfit_data.which_rows, false);
            
            % run optimization
            %
            results_ord = results_ord + 1;
            fprintf('\nFitting #%d with options:\n', results_ord);
            disp(options);
            
            results(results_ord) = mfit_optimize(likfun, param, mfit_data, nstarts);
            results_options(results_ord) = options;
            mfit_datas{results_ord} = mfit_data;

            fprintf('  BIC = %.4f, AIC = %.4f, params = \n', results(results_ord).bic, results(results_ord).aic);
            disp(results(results_ord).x);
        end
    end
end

if ~isempty(outfile)
    save(outfile, 'results', 'results_options');
end

toc
