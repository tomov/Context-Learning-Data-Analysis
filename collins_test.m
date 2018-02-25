function test_results = collins_test(stimuli, contexts, train_results, params, DO_PRINT)

% Collins & Frank clustering model. Clusters stimuli and contexts independently using DP (CRP).
%

Q0 = 0.1; % prior outcome expectaiton TODO collins = 0.5; Sam = 0; maybe fit?

eta = params(1); % learning rate
inv_softmax_temp = params(2); 
alpha = params(3); % concentration parameter
%alpha = 1.5; % TODO fix
 
if DO_PRINT
    disp('Collins Test');
    fprintf('eta = %.3f\n', eta);
    fprintf('beta = %.3f\n', inv_softmax_temp);
    fprintf('alpha = %.3f\n', alpha);
end

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

max_clusters = 5;

D_C = 3;          % # of contexts = K in Sam's model
D_S = size(stimuli, 2); % # of stimuli = D in Sam's model

N = size(stimuli, 1); % # observations

K_C = train_results.K_C; % # of active context clusters
K_S = train_results.K_S; % # of active stimulus clusters
P_C = train_results.P_C; % P_C(C,Z) = P(Z|C) cluster assignments for each context. Note row = C, col = Z
P_S = train_results.P_S; % P_S(S,Z) = P(Z|S) cluster assignments for each stimulus. Note row = S, col = Z

Q = train_results.Q; % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    s = find(stimuli(n,:));
    c = contexts(n); % one-hot vector for the context

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d  ============\n', n, c, s); end

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
    priors_C(:,:,n) = P_C;
    priors_S(:,:,n) = P_S;
    prior_c(n,:) = P_C(c,:);
    prior_s(n,:) = P_S(s,:);
    priors_Q(:,:,n) = Q;

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
end


test_results.choices = choices;
test_results.values = values;
test_results.priors_C = priors_C;
test_results.priors_S = priors_S;
test_results.prior_c = prior_c;
test_results.prior_s = prior_s;
test_results.priors_Q = priors_Q;

end

