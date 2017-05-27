function [choices, P_n, ww_n, P, ww, values, valuess, likelihoods, new_values, new_valuess, Sigma] = model_train(x, k, r, params, which_structures, DO_PRINT)

% Kalman filter to learn the context-cue-reward associations & posteriors
% for each causal structure
% See "Context-dependent learning and causal structure", Gershman, S.J. (2017)
%
% INPUT:
% x = matrix where each row is a stimulus vector for the given trial
% k = vector where each element is the context index on the given trial
% r = vector where each element is the outcome on the given trial
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

assert(numel(params) == 2);
prior_variance = params(1);
inv_softmax_temp = params(2);

momchil_mode = false; % include bias term for context in M2, so M2 = additive + modulatory hybrid

predict = @(V_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * V_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

% constants
%
N = size(x, 1); % # of trials
D = size(x, 2); % # of stimuli
K = 3;          % # of contexts

sigma_r = sqrt(0.01);
sigma_w = sqrt(prior_variance); % std for gaussian prior over weights, uncertainty; decreasing the learning rate
tau = sqrt(0.001);

% initialize Kalman filter
%
ww_n{1} = zeros(D, 1); % M1 weights: one per stimulus
ww_n{2} = zeros(D + 1, K); % M2 weights: one per stimulus-context pair
ww_n{3} = zeros(D + K, 1); % M3 weights: one per stimulus + one per context
ww_n{4} = zeros(D, 1); % M4 weights: one per context

Sigma_n{1} = sigma_w^2 * eye(D);
Sigma_n{2} = repmat(sigma_w^2 * eye(D + 1), 1, 1, K); % note the third dimension is the context
Sigma_n{3} = sigma_w^2 * eye(D + K);
Sigma_n{4} = sigma_w^2 * eye(K);

P_n = which_structures; % prior = 1 if we're including the model, 0 if not
P_n = P_n / sum(P_n);

% Store history for plotting and analysis
%
P = []; % history of posterior P(M | h_1:n)
ww{1} = []; % history of ww_1:n for M1
ww{2} = []; % history of ww_1:n for M2
ww{3} = []; % history of ww_1:n for M3
ww{4} = []; % history of ww_1:n for M4
choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)
valuess = []; % history of predicted outcomes, one for each model (causal structure)
likelihoods = []; % history of likelihoods, one for each model (causal structure)
new_values = []; % same as values but after the update (for the same stimulus)
new_valuess = []; % same as valuess but after the update (for the same stimulus)
Sigma{1} = []; % history of Sigma_1:n for M1
Sigma{2} = []; % history of Sigma_1:n for M2
Sigma{3} = []; % history of Sigma_1:n for M3
Sigma{4} = []; % history of Sigma_1:n for M4

% train
%
for n = 1:N % for each trial
    x_n = x(n, :)'; % stimulus at trial n
    k_n = k(n); % context idx at trial n
    r_n = r(n, :); % reward at trial n
    c_n = zeros(K, 1);
    c_n(k_n) = 1; % context vector like x_n
    if momchil_mode
        xb_n = [x_n; 1]; % stim vec + bias term
    else
        xb_n = [x_n; 0]; % stim vec + bias term
    end
    xx_n = [x_n; c_n]; % augmented stimulus + context vector
    
    % make a prediction based on h_1:n-1
    %
    value = @(x_n, xx_n, xb_n, k_n, c_n, ww_n, P_n) (x_n' * ww_n{1}) * P_n(1) + ... % M1 
                                         (xb_n' * ww_n{2}(:, k_n)) * P_n(2) + ... % M2
                                         (xx_n' * ww_n{3}) * P_n(3) + ... % M3   
                                         (c_n' * ww_n{4}) * P_n(4); % M4
    V_n = value(x_n, xx_n, xb_n, k_n, c_n, ww_n, P_n);
    out = predict(V_n);
    choices = [choices; out];
    values = [values; V_n];
    valuess = [valuess; x_n' * ww_n{1}, xb_n' * ww_n{2}(:, k_n), xx_n' * ww_n{3}, c_n' * ww_n{4}];
    
    % save(['kalman_state_', num2str(n), '.mat']);
    
    if DO_PRINT, fprintf('\npredction for x = %d, c = %d is %f (actual is %f)\n\n', find(x_n), c_n, out, r_n); end

    % get reward and update state
    %
    SSigma_n{1} = Sigma_n{1} + tau^2 * eye(D); % 1 / certainty for prediction x * w in M1
    SSigma_n{2} = Sigma_n{2}(:,:,k_n) + tau^2 * eye(D + 1); % 1 / certainty for prediction x * w in M2
    SSigma_n{3} = Sigma_n{3} + tau^2 * eye(D + K); % 1 / certainty for prediction x * w in M3
    SSigma_n{4} = Sigma_n{1} + tau^2 * eye(K); % 1 / certainty for prediction c * w in M4

    gain = @(x_n, SSigma_n) SSigma_n * x_n / (x_n' * SSigma_n * x_n + sigma_r^2);
    g_n{1} = gain(x_n, SSigma_n{1});    
    g_n{2} = gain(xb_n, SSigma_n{2});    
    g_n{3} = gain(xx_n, SSigma_n{3}); 
    g_n{4} = gain(c_n, SSigma_n{4});

    if DO_PRINT, fprintf('    g_ns = %.4f %.4f %.4f | %.4f %4.f %.4f | %.4f %.4f %.4f %.4f %.4f %.4f\n', g_n{1}, g_n{2}, g_n{3}); end

    Sigma_n{1} = SSigma_n{1} - g_n{1} * x_n' * SSigma_n{1};
    Sigma_n{2}(:,:,k_n) = SSigma_n{2} - g_n{2} * xb_n' * SSigma_n{2};
    Sigma_n{3} = SSigma_n{3} - g_n{3} * xx_n' * SSigma_n{3};
    Sigma_n{4} = SSigma_n{4} - g_n{4} * c_n' * SSigma_n{4};

    ww_n{1} = ww_n{1} + g_n{1} * (r_n - ww_n{1}' * x_n);
    ww_n{2}(:,k_n) = ww_n{2}(:,k_n) + g_n{2} * (r_n - ww_n{2}(:,k_n)' * xb_n);
    ww_n{3} = ww_n{3} + g_n{3} * (r_n - ww_n{3}' * xx_n);
    ww_n{4} = ww_n{4} + g_n{4} * (r_n - ww_n{4}' * c_n);

    if DO_PRINT
        disp('    ww_n{1} =');
        disp(ww_n{1});
        disp('    ww_n{2} =');
        disp(ww_n{2});
        disp('    ww_n{3} =');
        disp(ww_n{3});
        disp('    ww_n{4} =');
        disp(ww_n{4});
    end
    
    % likelihoods P(r_n | x_n, c_n, h_1:n-1, M)
    %
    likelihood = [normpdf(r_n, x_n' * ww_n{1}, x_n' * SSigma_n{1} * x_n + sigma_r^2), ...
                  normpdf(r_n, xb_n' * ww_n{2}(:,k_n), xb_n' * SSigma_n{2} * xb_n + sigma_r^2), ...
                  normpdf(r_n, xx_n' * ww_n{3}, xx_n' * SSigma_n{3} * xx_n + sigma_r^2), ...
                  normpdf(r_n, c_n' * ww_n{4}, c_n' * SSigma_n{4} * c_n + sigma_r^2)];
    likelihoods = [likelihoods; likelihood];
    
    % posteriors P(M | h_1:n)
    % 
    P_n(1) = P_n(1) * likelihood(1);
    P_n(2) = P_n(2) * likelihood(2);
    P_n(3) = P_n(3) * likelihood(3);
    P_n(4) = P_n(4) * likelihood(4);
    P_n = P_n / sum(P_n);
    if isnan(P_n(1))
        % edge case -- all likelihoods are 0, or the likelihoods of all
        % non-zero P_n's are 0 => the P_n's become all 0 after this => then
        % we divide by 0 and they all become NaNs
        % in to fix, interpret this scenario as them being all "equally
        % bad" => revert to the uniform
        %
        % TODO ask sam what to do
        %
        P_n = which_structures / sum(which_structures);
    end
    
    shit = [x_n' * ww_n{1}, xb_n' * ww_n{2}(:,k_n), xx_n' * ww_n{3}, c_n' * ww_n{4}];
    wtf = [x_n' * SSigma_n{1} * x_n + sigma_r^2, xb_n' * SSigma_n{2} * xb_n + sigma_r^2, xx_n' * SSigma_n{3} * xx_n + sigma_r^2, c_n' * SSigma_n{4} * c_n + sigma_r^2];
    pdfs = [normpdf(r_n, x_n' * ww_n{1}, x_n' * SSigma_n{1} * x_n + sigma_r^2), ...
        normpdf(r_n, xb_n' * ww_n{2}(:,k_n), xb_n' * SSigma_n{2} * xb_n + sigma_r^2), ...
        normpdf(r_n, xx_n' * ww_n{3}, xx_n' * SSigma_n{3} * xx_n + sigma_r^2), ...
        normpdf(r_n, c_n' * ww_n{4}, c_n' * SSigma_n{4} * c_n + sigma_r^2)];
    if DO_PRINT
        fprintf('trial %d\n', n);
        fprintf('       mean = %.4f %.4f %.4f %.4f\n', shit);
        fprintf('       var  = %.4f %.4f %.4f %.4f\n', wtf);
        fprintf('       pdfs = %.4f %.4f %.4f %.4f\n', pdfs);
        fprintf('       P    = %.4f %.4f %.4f %.4f\n', P_n);
    end

    P = [P; P_n];
    ww{1} = [ww{1}; ww_n{1}(1:2)'];
    ww{2} = [ww{2}; reshape(ww_n{2}(1:3,1:2), [1 6])];
    ww{3} = [ww{3}; ww_n{3}([1:2 4:5])'];
    ww{4} = [ww{4}; ww_n{4}(1:2)'];
    
    Sigma{1} = [Sigma{1}; Sigma_n{1}(eye(3) == 1)'];
    Sigma{2} = [Sigma{2}; Sigma_n{2}(eye(4) == 1)'];
    Sigma{3} = [Sigma{3}; [Sigma_n{3}(eye(6) == 1)', reshape(Sigma_n{3}([1:2 4:5], [1:2 4:5]), 16, 1)']];
    Sigma{4} = [Sigma{4}; Sigma_n{4}(eye(3) == 1)'];

    new_values = [new_values; value(x_n, xx_n, xb_n, k_n, c_n, ww_n, P_n)];
    new_valuess = [new_valuess; x_n' * ww_n{1}, xb_n' * ww_n{2}(:, k_n), xx_n' * ww_n{3}, c_n' * ww_n{4}];
    
 %   fprintf('            new Ps = %f %f %f\n', P(1), P(2), P(3));
end

assert(mean(sum(valuess(2:end,:) .* P(1:end-1,:), 2) == values(2:end)) == 1);
