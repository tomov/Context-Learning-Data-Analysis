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

contexts_onehot = zeros(N, D_C);
contexts_onehot(sub2ind([N D_C], 1:N, contexts')) = 1;
contexts = contexts_onehot;

K_C = train_results.K_C; % # of active context clusters
K_S = train_results.K_S; % # of active stimulus clusters
P_Zc_given_C = train_results.P_C; % P_C(C,Z) = P(Z|C) cluster assignments for each context. Note row = C, col = Z
P_Zs_given_S = train_results.P_S; % P_S(S,Z) = P(Z|S) cluster assignments for each stimulus. Note row = S, col = Z

Q = train_results.Q; % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    s = find(stimuli(n,:));
    c = find(contexts(n,:));

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d  ============\n', n, c, s); end

    % Initialize cluster probabilities according to CRP if stimulus and/or context is new
    %
    [P_Zc_given_C, K_C] = CRP_update(P_Zc_given_C, K_C, c, alpha, DO_PRINT);
    [P_Zs_given_S, K_S] = CRP_update(P_Zs_given_S, K_S, s, alpha, DO_PRINT);

    if DO_PRINT
        disp('      prior P(Z|C):');
        disp(P_Zc_given_C);
        disp('      prior P(Z|S):');
        disp(P_Zs_given_S);
        disp('      prior Q:');
        disp(Q);
    end
    priors_Zc_given_C(:,:,n) = P_Zc_given_C;
    priors_Zs_given_S(:,:,n) = P_Zs_given_S;
    prior_Zc_given_c(n,:) = P_Zc_given_C(:,c);
    prior_Zs_given_s(n,:) = P_Zs_given_S(:,s);
    priors_Q(:,:,n) = Q;

    % pick clusters of current stimulus/context for action selection (maximum a priori)
    %
    [~, z_c] = max(P_Zc_given_C(:,c));
    [~, z_s] = max(P_Zs_given_S(:,s));

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
test_results.priors_Zc_given_C = priors_Zc_given_C;
test_results.priors_Zs_given_S = priors_Zs_given_S;
test_results.prior_Zc_given_c = prior_Zc_given_c;
test_results.prior_Zs_given_s = prior_Zs_given_s;
test_results.priors_Q = priors_Q;

end


function [P, K] = CRP_update(P, K, c, alpha, DO_PRINT)
    % if context c is new, update its distribution over clusters P(Z_c|c)
    %
    if sum(P(:,c)) == 0
        % new context
        %
        p = [sum(P(1:K,:), 2); alpha];
        p = p / sum(p);
        assert(numel(p) == K+1);
        P(1:numel(p), c) = p;
        K = K + 1; % as per collins & frank's 2013 psych review paper
        if DO_PRINT, disp('                       New cluster'); end
        %if DO_PRINT, fprintf('      New : P(z|s) = [%s]\n', sprintf('%.3f, ', P_S(s,:))); end
    end
end

