function simulated = simulate_subjects(data, metadata, params, which_structures, which_rows, DO_PRINT)

% Simulate the subjects using the causal structure learning model, 
% given their stimulus sequences
% 
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% params = vector of hyperparameters:
%    params(:, 1) = prior_variance
%    params(:, 2) = inv_softmax_temp
%          for fixed effects, only have 1 row and these parameters will be
%          used for all subjects. For random effects, have a separate row
%          with the parameters for each subject
% which_structures = which causal structures to use, as a logical vector,
%                    e.g. [1 1 1 0] = M1, M2, M3
% which_rows (optional) = binary mask saying which rows from the data to
%                         use. If not supplied, assumes all rows are used.
% DO_PRINT (optional) = print output
%

if nargin < 5 || isempty(which_rows)
    which_rows = data.which_rows; % use all rows by default
else
    %warning('simulate_subjects: including rows that potentially contain "bad" subjects (%d total)', sum(which_rows));
    % never mind, we cannot include not good subjects -- we're iterating
    % over good subjects only b/c we've only fit parameters to them
    %
    assert(isequal(which_rows | data.which_rows, data.which_rows)); % make sure only good subjects are included
end
if nargin < 6 || isempty(DO_PRINT)
    DO_PRINT = true;
end

num_structures = 5;

% backwards compatibility: correct for extra structures from old fits
if numel(which_structures) < num_structures
    which_structures = [which_structures zeros(1, num_structures - numel(which_structures))];
    warning('using which_structures with less elements than there are structures');
end

simulated.keys = {}; % equivalent to response.keys but for the model (i.e. the responses)
simulated.pred = []; % the choice probability (not the actual choice) for each trial
simulated.P1 = []; % posterior P(M1 | ...) at each trial
simulated.P2 = []; % posterior P(M2 | ...) at each trial
simulated.P3 = []; % posterior P(M3 | ...) at each trial
simulated.P4 = []; % posterior P(M4 | ...) at each trial
simulated.P = []; % concatenated posteriors
simulated.Q1 = []; % prior Q(M1 | ...) at each trial
simulated.Q2 = []; % prior Q(M2 | ...) at each trial
simulated.Q3 = []; % prior Q(M3 | ...) at each trial
simulated.Q4 = []; % prior Q(M4 | ...) at each trial
simulated.Q = []; % concatenated priors
simulated.ww_after = {}; % weights for M1..M4 after update
simulated.ww_before = {}; % weights for M1..M4 before update 
simulated.ww_posterior = []; % weights after update for all structures concatenated
simulated.ww_prior = []; % weights before update for all structures concatenated
simulated.values = []; % values at each trial
simulated.valuess = []; % values for each sturcture at each trial
simulated.surprise = []; % D_KL at each trial for the posterior
simulated.KL_weights = []; % D_KL at each trial for the weights
simulated.likelihoods = []; % likelihoods for each causal structure
simulated.new_values = []; % updated values AFTER the trial
simulated.new_valuess = []; % updated values for each sturcture AFTER each trial
simulated.Sigma_after = {}; % Sigma for M1.M4 at each trial after update
simulated.Sigma_before = {}; % Sigma for M1..M4 at each trial before update
simulated.Sigma_posterior = []; % Sigma after update for all structures concatenated
simulated.Sigma_prior = []; % Sigmadefore update for all structures concatenated
simulated.lambdas = []; % lambdas at each trial

simulated.Q_values = []; % TODO hack fixme rm for q learning

% we either have the same params for all subjects (fixed effects)
% or a set of parameters for each subject (random effects)
%
assert(size(params, 1) == 1 || size(params, 1) == numel(metadata.subjects));
fixed_effects = size(params, 1) == 1;

s_id = 0;
for who = metadata.subjects
    s_id = s_id + 1;
    
    % figure out the parameters
    %
    if fixed_effects
        subject_params = params;
    else
        subject_params = params(s_id, :);
    end
    if DO_PRINT, fprintf('Simulating subj %s with params [%s]\n', who{1}, sprintf('%f ', subject_params)); end
    
    % Simulate each run separately
    %
    for condition = unique(data.contextRole)'
        which_runs = which_rows & strcmp(data.participant, who) & strcmp(data.contextRole, condition);
        runs = unique(data.runId(which_runs))';
        
        for run = runs
            which_train = which_runs & data.isTrain & data.runId == run;
            which_test = which_runs & ~data.isTrain & data.runId == run;

            % Get the training & test stimuli sequences for that run
            %
            [train_x, train_k, train_r, test_x, test_k] = convert_run(data, metadata, who, run);


            % TODO FIXME hack
            %
            if isequal(which_structures, 'simple_Q') % TODO fixme for q learning by reviewer 1

                train_results = q_train(train_x, train_k, train_r, subject_params, false);

                model_choices = train_results.choices > rand;
                model_response_keys = {};
                model_response_keys(model_choices) = {'left'};
                model_response_keys(~model_choices) = {'right'};
                simulated.keys(which_train) = model_response_keys;
                simulated.pred(which_train) = train_results.choices;
                simulated.values(which_train, :) = train_results.values;
                simulated.Q_values(:, :, which_train) = train_results.Qs;
           
                test_results = q_test(test_x, test_k, train_results, subject_params, false);

                model_test_choices = test_results.choices > rand;
                model_test_response_keys = {};
                model_test_response_keys(model_test_choices) = {'left'};
                model_test_response_keys(~model_test_choices) = {'right'};
                simulated.keys(which_test) = model_test_response_keys;
                simulated.pred(which_test) = test_results.choices;
                simulated.values(which_test, :) = test_results.values;

            elseif isequal(which_structures, 'Q_learning') % TODO fixme for q learning by reviewer 1

                train_results = q2_train(train_x, train_k, train_r, subject_params, false);

                model_choices = train_results.choices > rand;
                model_response_keys = {};
                model_response_keys(model_choices) = {'left'};
                model_response_keys(~model_choices) = {'right'};
                simulated.keys(which_train) = model_response_keys;
                simulated.pred(which_train) = train_results.choices;
                simulated.values(which_train, :) = train_results.values;
                simulated.Q_values(:, :, which_train) = train_results.Qs;
           
                test_results = q2_test(test_x, test_k, train_results, subject_params, false);

                model_test_choices = test_results.choices > rand;
                model_test_response_keys = {};
                model_test_response_keys(model_test_choices) = {'left'};
                model_test_response_keys(~model_test_choices) = {'right'};
                simulated.keys(which_test) = model_test_response_keys;
                simulated.pred(which_test) = test_results.choices;
                simulated.values(which_test, :) = test_results.values;

            elseif isequal(which_structures, 'simple_collins') % TODO fixme  reviewer 2

                train_results = collins_train(train_x, train_k, train_r, subject_params, false);

                model_choices = train_results.choices > rand;
                model_response_keys = {};
                model_response_keys(model_choices) = {'left'};
                model_response_keys(~model_choices) = {'right'};
                simulated.keys(which_train) = model_response_keys;
                simulated.pred(which_train) = train_results.choices;
                simulated.values(which_train, :) = train_results.values;
                simulated.priors_Zc_given_C(:,:,which_train) = train_results.priors_Zc_given_C;
                simulated.priors_Zs_given_S(:,:,which_train) = train_results.priors_Zs_given_S;
                simulated.prior_Zc_given_c(which_train,:) = train_results.prior_Zc_given_c;
                simulated.prior_Zs_given_s(which_train,:) = train_results.prior_Zs_given_s;
                simulated.priors_Q(:,:,which_train) = train_results.priors_Q;
                simulated.posteriors_Zc_given_C(:,:,which_train) = train_results.posteriors_Zc_given_C;
                simulated.posteriors_Zs_given_S(:,:,which_train) = train_results.posteriors_Zs_given_S;
                simulated.posterior_Zc_given_c(which_train,:) = train_results.posterior_Zc_given_c;
                simulated.posterior_Zs_given_s(which_train,:) = train_results.posterior_Zs_given_s;
                simulated.posteriors_Q(:,:,which_train) = train_results.posteriors_Q;
                simulated.PEs(which_train) = train_results.PEs;

                surprise_Zc_given_c = KL_divergence(train_results.posterior_Zc_given_c, train_results.prior_Zc_given_c);
                simulated.surprise_Zc_given_c(which_train, :) = surprise_Zc_given_c;

                surprise_Zs_given_s = KL_divergence(train_results.posterior_Zs_given_s, train_results.prior_Zs_given_s);
                simulated.surprise_Zs_given_s(which_train, :) = surprise_Zs_given_s;

                test_results = collins_test(test_x, test_k, train_results, subject_params, false);

                model_test_choices = test_results.choices > rand;
                model_test_response_keys = {};
                model_test_response_keys(model_test_choices) = {'left'};
                model_test_response_keys(~model_test_choices) = {'right'};
                simulated.keys(which_test) = model_test_response_keys;
                simulated.pred(which_test) = test_results.choices;
                simulated.values(which_test, :) = test_results.values;
                simulated.priors_Zc_given_C(:,:,which_test) = test_results.priors_Zc_given_C;
                simulated.priors_Zs_given_S(:,:,which_test) = test_results.priors_Zs_given_S;
                simulated.prior_Zc_given_c(which_test,:) = test_results.prior_Zc_given_c;
                simulated.prior_Zs_given_s(which_test,:) = test_results.prior_Zs_given_s;
                simulated.priors_Q(:,:,which_test) = test_results.priors_Q;
                simulated.posteriors_Zc_given_C(:,:,which_test) = test_results.priors_Zc_given_C;
                simulated.posteriors_Zs_given_S(:,:,which_test) = test_results.priors_Zs_given_S;
                simulated.posterior_Zc_given_c(which_test,:) = test_results.prior_Zc_given_c;
                simulated.posterior_Zs_given_s(which_test,:) = test_results.prior_Zs_given_s;
                simulated.posteriors_Q(:,:,which_test) = test_results.priors_Q;

            else
                assert(~ischar(which_structures));
                % For a given run of a given subject, run the model on the same
                % sequence of stimuli and see what it does.
                %
                %[choices, P_n, ww_n, P, ww_after, values, valuess, likelihoods, new_values, new_valuess, Sigma_after, lambdas, ww_before, Sigma_before] = ...
                train_results = model_train(train_x, train_k, train_r, subject_params, which_structures, false);

                model_choices = train_results.choices > rand;
                model_response_keys = {};
                model_response_keys(model_choices) = {'left'};
                model_response_keys(~model_choices) = {'right'};
                simulated.keys(which_train) = model_response_keys;
                simulated.pred(which_train) = train_results.choices;
                simulated.P1(which_train) = train_results.P(:,1);
                simulated.P2(which_train) = train_results.P(:,2);
                simulated.P3(which_train) = train_results.P(:,3);
                simulated.P4(which_train) = train_results.P(:,4);
                simulated.P(which_train, :) = train_results.P;
                priors = which_structures / sum(which_structures);
                Q = [priors; train_results.P(1:end-1,:)];
                simulated.Q1(which_train) = Q(:,1);
                simulated.Q2(which_train) = Q(:,2);
                simulated.Q3(which_train) = Q(:,3);
                simulated.Q4(which_train) = Q(:,4);
                simulated.Q(which_train, :) = Q;
                simulated.ww_after{1}(which_train, :) = train_results.ww_after{1};
                simulated.ww_after{2}(which_train, :) = train_results.ww_after{2};
                simulated.ww_after{3}(which_train, :) = train_results.ww_after{3};
                simulated.ww_after{4}(which_train, :) = train_results.ww_after{4};
                simulated.ww_before{1}(which_train, :) = train_results.ww_before{1};
                simulated.ww_before{2}(which_train, :) = train_results.ww_before{2};
                simulated.ww_before{3}(which_train, :) = train_results.ww_before{3};
                simulated.ww_before{4}(which_train, :) = train_results.ww_before{4};
                simulated.values(which_train, :) = train_results.values;
                simulated.valuess(which_train, :) = train_results.valuess;
                surprise = KL_divergence(train_results.P, Q);
                simulated.surprise(which_train, :) = surprise;
                simulated.likelihoods(which_train, :) = train_results.likelihoods;
                simulated.new_values(which_train, :) = train_results.new_values;
                simulated.new_valuess(which_train, :) = train_results.new_valuess;
                simulated.Sigma_after{1}(:, :, which_train) = train_results.Sigma_after{1};
                simulated.Sigma_after{2}(:, :, which_train) = train_results.Sigma_after{2};
                simulated.Sigma_after{3}(:, :, which_train) = train_results.Sigma_after{3};
                simulated.Sigma_before{1}(:, :, which_train) = train_results.Sigma_before{1};
                simulated.Sigma_before{2}(:, :, which_train) = train_results.Sigma_before{2};
                simulated.Sigma_before{3}(:, :, which_train) = train_results.Sigma_before{3};
                simulated.lambdas(which_train, :) = train_results.lambdas;
                
                % Concatenate weights 
                %
                dim = 0;
                for M = 1:3
                    dim = dim + size(train_results.ww_before{M}, 2); 
                end
                simulated.ww_prior(which_train, :) = zeros(sum(which_train), dim);
                simulated.ww_posterior(which_train, :) = zeros(sum(which_train), dim);
                simulated.Sigma_prior(:, :, which_train) = zeros(dim, dim, sum(which_train));
                simulated.Sigma_posterior(:, :, which_train) = zeros(dim, dim, sum(which_train));
                
                dim = 0;
                for M = 1:3
                    dims = dim + 1 : dim + size(train_results.ww_before{M}, 2);
                    
                    simulated.ww_prior(which_train, dims) = train_results.ww_before{M};
                    simulated.Sigma_prior(dims, dims, which_train) = train_results.Sigma_before{M};
                    
                    simulated.ww_posterior(which_train, dims) = train_results.ww_after{M};
                    simulated.Sigma_posterior(dims, dims, which_train) = train_results.Sigma_after{M};
                    
                    dim = dim + size(train_results.ww_before{M}, 2); 
                end
                

                % See what the model predicts for the test trials of that run
                %
                %[test_choices, test_values, test_valuess] = ...
                test_results = model_test(test_x, test_k, train_results, subject_params);

                model_test_choices = test_results.choices > rand;
                model_test_response_keys = {};
                model_test_response_keys(model_test_choices) = {'left'};
                model_test_response_keys(~model_test_choices) = {'right'};
                simulated.keys(which_test) = model_test_response_keys;
                simulated.pred(which_test) = test_results.choices;
                simulated.values(which_test, :) = test_results.values;
                simulated.valuess(which_test, :) = test_results.valuess;

                % Get the subject's responses too.
                %
                resp = data.response.keys(which_train);
                human_choices = strcmp(resp, 'left'); % sick == 1            

                % compute the KL divergence for the weights
                %
                ww_prior = [train_results.ww_before{1} train_results.ww_before{2} train_results.ww_before{3}];
                Sigma_prior = nan(size(ww_prior, 2), size(ww_prior, 2), size(ww_prior, 1));
                for i = 1:size(Sigma_prior, 3)
                    Sigma_prior(:,:,i) = blkdiag(train_results.Sigma_before{1}(:,:,i), train_results.Sigma_before{2}(:,:,i), train_results.Sigma_before{3}(:,:,i));
                end
                ww_posterior = [train_results.ww_after{1} train_results.ww_after{2} train_results.ww_after{3}];
                Sigma_posterior = nan(size(ww_posterior, 2), size(ww_posterior, 2), size(ww_posterior, 1));
                for i = 1:size(Sigma_posterior, 3)
                    Sigma_posterior(:,:,i) = blkdiag(train_results.Sigma_after{1}(:,:,i), train_results.Sigma_after{2}(:,:,i), train_results.Sigma_after{3}(:,:,i));
                end
                simulated.KL_weights(which_train, :) = KL_divergence_gauss(ww_posterior, Sigma_posterior, ww_prior, Sigma_prior);

            end % end TODO fixme hack

        end
    end
end

simulated.keys = simulated.keys';
simulated.pred = simulated.pred';
simulated.P1 = simulated.P1';
simulated.P2 = simulated.P2';
simulated.P3 = simulated.P3';
simulated.P4 = simulated.P4';
simulated.Q1 = simulated.Q1';
simulated.Q2 = simulated.Q2';
simulated.Q3 = simulated.Q3';
simulated.Q4 = simulated.Q4';




