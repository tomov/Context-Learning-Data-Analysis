function test_results = q2_test(x, k, train_results, params, DO_PRINT)
% same as model_test but for Q-learning model proposed by reviewer #1, except w learning rate

% INPUT:
% x = matrix where each row is a stimulus vector for the given trial
% k = vector where each element is the context index on the given trial
% 

assert(numel(params) == 2 || numel(params) == 3);
learning_rate = params(1);
inv_softmax_temp = params(2);

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

% constants
%
N = size(x, 1); % # of trials
D = size(x, 2); % # of stimuli
K = 3;          % # of contexts


% Store history for plotting and analysis
%
choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

% test 
%
for n = 1:N % for each trial
    x_n = x(n, :)'; % stimulus at trial n
    k_n = k(n); % context idx at trial n

    s = k_n; % stimulus = context
    assert(numel(s) == 1);
    a = find(x_n); % action = cue

    % make a prediction based on h_1:n-1
    %
    if train_results.count(s,a) == 0
        % no prior experience with this cue/context
        %
        %assert(sum(count(s,:)) == 0 || sum(count(:,a) == 0));
        if sum(train_results.count(s,:)) == 0 && sum(train_results.count(:,a)) == 0 % no prior experience with cue nor context
            V_n = 0.5; % TODO prior prediciton -- free parameter
        else % prior experience with either cue or context (or both) -> average over them
            V_n = sum(train_results.Q(:,a) .* train_results.count(:,a)) + sum(train_results.Q(s,:) .* train_results.count(s,:));
            V_n = V_n / (sum(train_results.count(:,a)) + sum(train_results.count(s,:)));
        end
    else
        % seen it
        %
        V_n = train_results.Q(s,a);
    end

    out = predict(V_n);
    choices = [choices; out];
    values = [values; V_n];

    if DO_PRINT
        fprintf('trial %d: x(a) = %d, k(s) = %d;  V_n = %.2f, out = %.2f\n', n, a, k_n, V_n, out);
        disp(train_results.Q);
    end
end


test_results.choices = choices;
test_results.values = values;

