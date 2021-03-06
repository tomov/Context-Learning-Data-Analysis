function train_results = model_train(stimuli, contexts, rewards, params, which_structures, DO_PRINT)

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
w{1} = initial_weight * ones(D, 1); % M1 weights: one per stimulus
w{2} = initial_weight * ones(D, K); % M2 weights: one per stimulus-context pair
w{3} = initial_weight * ones(D + K, 1); % M3 weights: one per stimulus + one per context
w{4} = initial_weight * ones(K, 1); % M4 weights: one per context
w{5} = initial_weight * ones(K, D); % M4 weights: one per stimulus-context pair

S{1} = sigma_w^2 * eye(D);
S{2} = repmat(sigma_w^2 * eye(D), 1, 1, K); % note the third dimension is the context
S{3} = sigma_w^2 * eye(D + K);
S{4} = sigma_w^2 * eye(K);
S{5} = repmat(sigma_w^2 * eye(K), 1, 1, D); % note the third dimension is the cue

P = which_structures / sum(which_structures);
assert(numel(P) == num_structures);

% Store history for plotting and analysis
%
Posterior = []; % history of posterior P(M | h_1:n)
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


% train
%
for n = 1:N % for each trial
    % Make prediction for current trial
    % notice conversion to column vectors
    %
    r = rewards(n);
    c = full(ind2vec(contexts(n), K)); % one-hot vector for the context
    x{1} = stimuli(n,:)';
    x{2} = x{1};
    x{3} = [x{1}; c];
    x{4} = c;
    x{5} = c;
 
    [V, vals] = model_value_helper(x, w, stimuli(n,:), contexts(n), P);

    out = predict(V);
    choices = [choices; out];
    values = [values; V];
    valuess = [valuess; vals];
   
    % log 
    %
    Sigma_before{1}(:,:,n) = S{1}(1:2,1:2);
    Sigma_before{2}(1:2,1:2,n) = S{2}(1:2,1:2,1);
    Sigma_before{2}(3:4,3:4,n) = S{2}(1:2,1:2,2);
    Sigma_before{3}(:,:,n) = S{3}([1 2 4 5], [1 2 4 5]);
    Sigma_before{4}(:,:,n) = S{4}(1:2,1:2);
    ww_before{1} = [ww_before{1}; w{1}(1:2)'];
    ww_before{2} = [ww_before{2}; reshape(w{2}(1:2,1:2), [1 4])];
    ww_before{3} = [ww_before{3}; w{3}([1:2 4:5])'];
    ww_before{4} = [ww_before{4}; w{4}(1:2)'];

    % Kalman filter update 
    %
    for i = 1:num_structures
        if i == 2
            % for M2, each context has a separate Kalman filter
            % Notice we update all of them but only the active one 'sees' the real stimulus, and we get the prediction from it
            %
            for k = 1:K 
                if k == contexts(n)
                    H = x{i}';
                else
                    H = zeros(1, D); % inactive contexts updated with empty stimulus -> only covariance changes
                end
                F = eye(size(w{i}(:,k), 1));
                B = zeros(size(F));
                u = zeros(size(w{i}(:,k)));
                Q = eye(size(w{i}(:,k), 1)) * tau^2;
                [w{i}(:,k), S{i}(:,:,k), r_mean_temp, r_var_temp] = kalman(r, w{i}(:,k), S{i}(:,:,k), F, B, u, Q, H, sigma_r^2);
                if k == contexts(n) % only consider prediction for current context
                    r_mean{i} = r_mean_temp;
                    r_var{i} = r_var_temp;
                end
            end
        elseif i == 5
            % for M5 = M2', each cue has a separate Kalman filter
            % Notice we update all of them but only the active one 'sees' the real context, and we get the prediction from it
            %
            for d = 1:D
                if stimuli(n,d)
                    H = x{i}';
                else
                    H = zeros(1, K); % inactive cues updated with empty stimulus -> only covariance changes
                end
                F = eye(size(w{i}(:,d), 1));
                B = zeros(size(F));
                u = zeros(size(w{i}(:,d)));
                Q = eye(size(w{i}(:,d), 1)) * tau^2;
                [w{i}(:,d), S{i}(:,:,d), r_mean_temp, r_var_temp] = kalman(r, w{i}(:,d), S{i}(:,:,d), F, B, u, Q, H, sigma_r^2);
                if stimuli(n,d)
                    r_mean{i} = r_mean_temp;
                    r_var{i} = r_var_temp;
                end
            end
        else
            F = eye(size(w{i}, 1));
            B = zeros(size(F));
            u = zeros(size(w{i}));
            Q = eye(size(w{i}, 1)) * tau^2;
            [w{i}, S{i}, r_mean{i}, r_var{i}] = kalman(r, w{i}, S{i}, F, B, u, Q, x{i}', sigma_r^2);
        end 
    end

    % Update posterior over causal structures
    %
    for i = 1:num_structures
        liks(i) = normpdf(r, r_mean{i}, r_var{i});
    end
    P = P .* liks;
    P = P / sum(P);
    if isnan(P(1))
        % edge case -- all likelihoods are 0, or the likelihoods of all
        % non-zero P_n's are 0 => the P_n's become all 0 after this => then
        % we divide by 0 and they all become NaNs
        % in to fix, interpret this scenario as them being all "equally
        % bad" => revert to the uniform
        %
        % TODO SAM ask sam what to do
        %
        P = which_structures / sum(which_structures);
    end

    % log stuff
    %
    likelihoods = [likelihoods; liks];
    Posterior = [Posterior; P];
    lambdas = [lambdas; r_var{1} r_var{2} r_var{3} r_var{4} r_var{5}];
    ww_after{1} = [ww_after{1}; w{1}(1:2)'];
    ww_after{2} = [ww_after{2}; w{2}(1:2,1)' w{2}(1:2,2)'];
    ww_after{3} = [ww_after{3}; w{3}([1:2 4:5])'];
    ww_after{4} = [ww_after{4}; w{4}(1:2)'];
    Sigma_after{1}(:,:,n) = S{1}(1:2,1:2);
    Sigma_after{2}(1:2,1:2,n) = S{2}(1:2,1:2,1);
    Sigma_after{2}(3:4,3:4,n) = S{2}(1:2,1:2,2);
    Sigma_after{3}(:,:,n) = S{3}([1 2 4 5], [1 2 4 5]);
    Sigma_after{4}(:,:,n) = S{4}(1:2,1:2);
    [V, vals] = model_value_helper(x, w, stimuli(n,:), contexts(n), P);
    new_values = [new_values; V];
    new_valuess = [new_valuess; vals];
end

 
assert(mean(sum(valuess(2:end,:) .* Posterior(1:end-1,:), 2) - values(2:end)) < 1e-6);


train_results.choices = choices;
train_results.P_n = P;
train_results.ww_n = w;
train_results.P = Posterior;
train_results.ww_after = ww_after;
train_results.values = values;
train_results.valuess = valuess;
train_results.likelihoods = likelihoods;
train_results.new_values = new_values;
train_results.new_valuess = new_valuess;
train_results.Sigma_after = Sigma_after;
train_results.lambdas = lambdas; 
train_results.ww_before = ww_before;
train_results.Sigma_before = Sigma_before;

end

