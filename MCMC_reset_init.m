function particle=ideal_init(stimuli, contexts, rewards, params, which_structures, DO_PRINT)

    % copy of ideal_init

% Kalman filter to learn the context-cue-reward associations & posteriors
% for each causal structure
% See "Context-dependent learning and causal structure", Gershman, S.J. (2017)
%
% Causal structures:
% M1 = irrelevant context
% M2 = modulatory context
% M3 = additive context
% M4 = irrelevant cue = M1' of reviewer 2
% M5 = modulatory cue = M2' of reviewer 2
%
% INPUT:
% stimuli = x = matrix where each row is a stimulus vector for the given trial
% contexts = k = vector where each element is the context index on the given trial
% rewards = r = vector where each element is the outcome on the given trial
% params = vector of hyperparameters
%    params(1) = prior_variance
%    params(2) = inv_softmax_temp
% which_structures = logical vector indicating which causal structures to
%                    use in the model, e.g. [1 1 1 0] for M1, M2, M3
% DO_PRINT = whether to print stuff for debugging
% 
% OUTPUT:
% choices = choice probability on each trial
% P_n = posterior over causal structures after training
% ww_n = cell array with the final learned weights for each causal structure
% P = history of posteriors over causal structures (columns) after each
%     trial (rows)
% ww = history of ww_n after each trial
% values = the expected outcome V_n after each trial
% valuess = the expected outcome predicted by each causal structure
%           (columns) after each trial (rows)
% likelihoods = the likelihood of the outcome on each trial (rows) for each
%               causal structure (columns)
% new_values = expected outcome V_n as predicted after seeing the outcome 
%              (and updating the parameters). I.e. like values but after
%              the update
% new_valuess = expected outcome as predicted by each causal structure
%               (columns) after seeing the outcome on the given trial
%               (rows). I.e. like valuess but after the update
% Sigma = history of Sigma's (weight covariance matrices) for each causal
%         structure

%disp(params);

assert(numel(params) == 2 || numel(params) == 3 || numel(params) == 4);
prior_variance = params(1);
inv_softmax_temp = params(2);
if numel(params) >= 3
    diffusion_variance = params(3);
else
    diffusion_variance = 0.001;
end
if numel(params) >= 4
    initial_weight = params(4);
else
    initial_weight = 0;
end

predict = @(V_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * V_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

% constants
%
N = size(stimuli, 1); % # of trials
D = size(stimuli, 2); % # of stimuli
K = 3;          % # of contexts
num_structures = 5;

sigma_r = sqrt(0.01);
sigma_w = sqrt(prior_variance); % std for gaussian prior over weights, uncertainty; decreasing the learning rate
tau = sqrt(diffusion_variance);

% initialize Kalman filter
%
particle.w{1} = initial_weight * ones(D, 1); % M1 weights: one per stimulus
particle.w{2} = initial_weight * ones(D, K); % M2 weights: one per stimulus-context pair
particle.w{3} = initial_weight * ones(D + K, 1); % M3 weights: one per stimulus + one per context
particle.w{4} = initial_weight * ones(K, 1); % M4 weights: one per context
particle.w{5} = initial_weight * ones(K, D); % M4 weights: one per stimulus-context pair
particle.S{1} = sigma_w^2 * eye(D);
particle.S{2} = repmat(sigma_w^2 * eye(D), 1, 1, K); % note the third dimension is the context
particle.S{3} = sigma_w^2 * eye(D + K);
particle.S{4} = sigma_w^2 * eye(K);
particle.S{5} = repmat(sigma_w^2 * eye(K), 1, 1, D); % note the third dimension is the cue
particle.P = which_structures / sum(which_structures);

particle.sample = zeros(size(particle.P));
ix = randsample(length(particle.P), 1, true, particle.P);
particle.sample(ix) = 1;

particle.P_prop = which_structures / sum(which_structures);

particle.Posterior = [particle.P]; % history of posterior P(M | h_1:n)
particle.samples = [particle.sample]; % history of samples
%{
% Store history for plotting and analysis
%
ww_after{1} = []; % history of ww_1:n for M1
ww_after{2} = []; % history of ww_1:n for M2
ww_after{3} = []; % history of ww_1:n for M3
ww_after{4} = []; % history of ww_1:n for M4
ww_before{1} = []; % history of ww_1:n for M1 before the update
ww_before{2} = []; % history of ww_1:n for M2 before the update
ww_before{3} = []; % history of ww_1:n for M3 before the update
ww_before{4} = []; % history of ww_1:n for M4 before the update
choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)
valuess = []; % history of predicted outcomes, one for each model (causal structure)
likelihoods = []; % history of likelihoods, one for each model (causal structure)
new_values = []; % same as values but after the update (for the same stimulus)
new_valuess = []; % same as valuess but after the update (for the same stimulus)
Sigma_after{1} = zeros(2, 2, N); % history of Sigma_1:n for M1 after
Sigma_after{2} = zeros(4, 4, N); % history of Sigma_1:n for M2 after
Sigma_after{3} = zeros(4, 4, N); % history of Sigma_1:n for M3 after
Sigma_after{4} = []; % history of Sigma_1:n for M4 after
Sigma_before{1} = zeros(2, 2, N); % history of Sigma_1:n for M1 before
Sigma_before{2} = zeros(4, 4, N); % history of Sigma_1:n for M2 before
Sigma_before{3} = zeros(4, 4, N); % history of Sigma_1:n for M3 before
Sigma_before{4} = []; % history of Sigma_1:n for M4 before
lambdas = []; % history of lambdas
%}
