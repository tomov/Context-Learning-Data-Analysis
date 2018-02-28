function test_results = flat_test(stimuli, contexts, train_results, params, DO_PRINT)

% Collins & Frank flat RL model
%

eta = params(1); % learning rate
inv_softmax_temp = params(2); 
%alpha = 1.5; % TODO fix
 
if DO_PRINT
    disp('Collins Test');
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

Q = train_results.Q; % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    s = find(stimuli(n,:));
    c = find(contexts(n,:));

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d  ============\n', n, c, s); end

    % pick action
    %
    q = Q(c, s);
    out = predict(q);
    choices = [choices; out];
    values = [values; q];

    Qs(:,:,n) = Q;

    if DO_PRINT, fprintf('   choice = %.3f\n', out); end
end


test_results.choices = choices;
test_results.values = values;
test_results.Qs = Qs;
