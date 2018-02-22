function train_results = collins_train(stimuli, contexts, rewards, params, DO_PRINT)

% Collins & Frank clustering model. Clusters stimuli and contexts independently using DP (CRP).
%

Q0 = 0.5; % prior outcome expectaiton TODO collins = 0.5; Sam = 0; maybe fit?

eta = params(1); % learning rate
inv_softmax_temp = params(2); 
alpha = params(3); % concentration parameter

if DO_PRINT
    disp('Collins Train');
    fprintf('eta = %.3f\n', eta);
    fprintf('beta = %.3f\n', inv_softmax_temp);
    fprintf('alpha = %.3f\n', alpha);
end

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

max_clusters = 5;

D_C = 3;          % # of contexts = K in Sam's model
D_S = size(stimuli, 2); % # of stimuli = D in Sam's model

N = size(stimuli, 1); % # observations


K_C = 0; % # of active context clusters
K_S = 0; % # of active stimulus clusters
P_C = zeros(D_C, max_clusters); % P_C(C,Z) = P(Z|C) cluster assignments for each context. Note row = C, col = Z
P_S = zeros(D_S, max_clusters); % P_S(S,Z) = P(Z|S) cluster assignments for each stimulus. Note row = S, col = Z

Q = Q0 * ones(max_clusters, max_clusters); % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    r = rewards(n);
    s = find(stimuli(n,:));
    c = contexts(n); % one-hot vector for the context

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d -> r = %d ============\n', n, c, s, r); end

    % Initialize cluster probabilities according to CRP if stimulus and/or context is new
    %
    if sum(P_C(c,:)) == 0
        % new context
        %
        p = [sum(P_C(:,1:K_C)), alpha];
        p = p / sum(p);
        P_C(c,1:numel(p)) = p;
        K_C = K_C + 1; % as per their 2013 psych review paper
        if DO_PRINT, fprintf('      New c: P(z|c) = [%s]\n', sprintf('%.3f, ', P_C(c,:))); end
    end

    if sum(P_S(s,:)) == 0
        % new stimulus 
        %
        p = [sum(P_S(:,1:K_S)), alpha];
        p = p / sum(p);
        P_S(s,1:numel(p)) = p;
        K_S = K_S + 1;
        if DO_PRINT, fprintf('      New s: P(z|s) = [%s]\n', sprintf('%.3f, ', P_S(s,:))); end
    end

    if DO_PRINT
        disp('      prior P(Z|C):');
        disp(P_C);
        disp('      prior P(Z|S):');
        disp(P_S);
        disp('      prior Q:');
        disp(Q);
    end


    % pick clusters of current stimulus/context for action selection (maximum a priori)
    %
    [~, z_c] = max(P_C(c,:));
    [~, z_s] = max(P_S(s,:));

    if DO_PRINT, fprintf('   (a priori) inferred z_c = %d, z_s = %d, Q = %.4f\n', z_c, z_s, Q(z_c, z_s)); end

    % pick action
    %
    q = Q(z_c, z_s);
    out = predict(q);
    choices = [choices; out];
    values = [values; q];

    if DO_PRINT, fprintf('   choice = %.3f\n', out); end

    % update posterior over clusters for current stimulus and context
    %
    assert(r == 1 || r == 0); % assumes binary rewards

    if DO_PRINT, fprintf('          across context clusters: Q(:, z_s) = %s\n', sprintf('%.3f, ', Q(:,z_s))); end
    for i=1:K_C % for each context cluster
        q = Q(i, z_s);
        lik = binopdf(r, 1, q);
        if DO_PRINT, fprintf('                 cluster %d: Q(%d, z_s) = %.3f, lik = %.3f, prior P(%d|c) = %.3f\n', i, i, q, lik, i, P_C(c, i)); end
        P_C(c, i) = P_C(c, i) * lik;
    end
    P_C(c, :) = P_C(c, :) / sum(P_C(c, :));
    if DO_PRINT, fprintf('      Updated c: P(z|c) = [%s]\n', sprintf('%.3f, ', P_C(c,:))); end

    if DO_PRINT, fprintf('          across stimulus clusters: Q(z_c, :) = %s\n', sprintf('%.3f, ', Q(z_c,:))); end
    for i=1:K_S % for each stimulus cluster
        q = Q(z_c, i);
        lik = binopdf(r, 1, q);
        if DO_PRINT, fprintf('                 cluster %d: Q(z_c, %d) = %.3f, lik = %.3f, prior P(%d|s) = %.3f\n', i, i, q, lik, i, P_S(s, i)); end
        P_S(s, i) = P_S(s, i) * lik;
    end
    P_S(s, :) = P_S(s, :) / sum(P_S(s, :));
    if DO_PRINT, fprintf('      Updated s: P(z|s) = [%s]\n', sprintf('%.3f, ', P_S(s,:))); end

    % pick clusters of current stimulus/context for updating (maximum a posteriori)
    %
    [~, z_c] = max(P_C(c,:));
    [~, z_s] = max(P_S(s,:));

    if DO_PRINT, fprintf('   (a posteriori) inferred z_c = %d, z_s = %d, Q = %.4f\n', z_c, z_s, Q(z_c, z_s)); end

    % update Q-value based on reward
    %
    PE = r - Q(z_c, z_s);
    Q(z_c, z_s) = Q(z_c, z_s) + eta * PE;

    if DO_PRINT, fprintf('   PE = %.3f - %.3f = %.3f\n', r, Q(z_c, z_s), PE); end
    if DO_PRINT, fprintf('     new Q(%d, %d) = %.3f\n', z_c, z_s, Q(z_c, z_s)); end

    if DO_PRINT
        disp('      posterior P(Z|C):');
        disp(P_C);
        disp('      posterior P(Z|S):');
        disp(P_S);
        disp('      posterior Q:');
        disp(Q);
    end

end


train_results.choices = choices;
train_results.values = values;
train_results.P_C = P_C;
train_results.P_S = P_S;
train_results.K_C = K_C;
train_results.K_S = K_S;
train_results.Q = Q;
