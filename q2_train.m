function train_results = q2_train(x, k, r, params, DO_PRINT)
% same as model_train but for Q-learning model proposed by reviewer #1, + learning rate

% INPUT:
% x = matrix where each row is a stimulus vector for the given trial
% k = vector where each element is the context index on the given trial
% r = vector where each element is the outcome on the given trial
% 

assert(numel(params) == 2);
learning_rate = params(1);
inv_softmax_temp = params(2);

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

% constants
%
N = size(x, 1); % # of trials
D = size(x, 2); % # of stimuli
K = 3;          % # of contexts

alpha = learning_rate;

% initialize Q learning
%
Q = zeros(K,D); % analogous to our initialization of the kalman filter weights to 0, i.e. predict no sickness by default
count = zeros(K,D);

% Store history for plotting and analysis
%

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)
Qs = nan(K,D,N); % history of Q-values

% train
%
for n = 1:N % for each trial
    x_n = x(n, :)'; % stimulus at trial n
    k_n = k(n); % context idx at trial n
    r_n = r(n, :); % reward at trial n

    s = k_n; % stimulus = context
    assert(numel(s) == 1);
    a = find(x_n); % action = cue

    % make a prediction based on h_1:n-1
    %
    if count(s,a) == 0
        % no prior experience with this cue/context
        %
        %assert(sum(count(s,:)) == 0 || sum(count(:,a) == 0));
        if sum(count(s,:)) == 0 && sum(count(:,a)) == 0 % no prior experience with cue nor context
            V_n = 0; % predict no sickness for unseen anything (only happens on first trial; analogous to our zero initialized kalman weights)
        else % prior experience with either cue or context (or both) -> average over them
            V_n = sum(Q(:,a) .* count(:,a)) + sum(Q(s,:) .* count(s,:));
            V_n = V_n / (sum(count(:,a)) + sum(count(s,:)));
        end
    else
        % seen it
        %
        V_n = Q(s,a);
    end

    out = predict(V_n);
    choices = [choices; out];
    values = [values; V_n];

    if DO_PRINT
        fprintf('trial %d: x(a) = %d, k(s) = %d, r = %d;  V_n = %.2f, out = %.2f\n', n, a, k_n, r_n, V_n, out);
        disp(Q);
    end

    % get reward and update state
    %
    count(s,a) = count(s,a) + 1;

    Q(s,a) = Q(s,a) + alpha * (r_n - Q(s,a));

    Qs(:,:,n) = Q; 
end



train_results.choices = choices;
train_results.values = values;
train_results.Qs = Qs;
train_results.Q = Q;
train_results.count = count;
