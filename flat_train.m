function train_results = flat_train(stimuli, contexts, rewards, params, DO_PRINT)

% Collins & Frank flat RL model.
%

assert(numel(params) == 2 || numel(params) == 3);
eta = params(1); % learning rate
inv_softmax_temp = params(2); 
if numel(params) >= 3
    Q0 = params(3);
else
    Q0 = 0.1; % prior outcome expectaiton; collins = 0.5; Sam = 0;
end

if DO_PRINT
    disp('Collins Train');
    fprintf('eta = %.3f\n', eta);
    fprintf('beta = %.3f\n', inv_softmax_temp);
end

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

max_clusters = 5;

D_C = 3;          % # of contexts = K in Sam's model
D_S = size(stimuli, 2); % # of stimuli = D in Sam's model

N = size(stimuli, 1); % # observations

contexts_onehot = zeros(N, D_C);
contexts_onehot(sub2ind([N D_C], 1:N, contexts')) = 1;
contexts = contexts_onehot;


Q = Q0 * ones(max_clusters, max_clusters); % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    r = rewards(n);
    s = find(stimuli(n,:));
    c = find(contexts(n,:));

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d -> r = %d ============\n', n, c, s, r); end

    % pick action
    %
    q = Q(c, s);
    out = predict(q);
    choices = [choices; out];
    values = [values; q];

    if DO_PRINT, fprintf('   choice = %.3f\n', out); end

    % update Q-value based on reward
    %
    PE = r - Q(c, s);
    Q(c, s) = Q(c, s) + eta * PE;

    PEs(n) = PE;
    Qs(:,:,n) = Q;
end


train_results.choices = choices;
train_results.values = values;
train_results.Q = Q;
train_results.Qs = Qs;
train_results.PEs = PEs;
