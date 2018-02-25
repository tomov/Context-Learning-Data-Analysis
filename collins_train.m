function train_results = collins_train(stimuli, contexts, rewards, params, DO_PRINT)

% Collins & Frank clustering model. Clusters stimuli and contexts independently using DP (CRP).
%

Q0 = 0.1; % prior outcome expectaiton TODO collins = 0.5; Sam = 0; maybe fit?

eta = params(1); % learning rate
inv_softmax_temp = params(2); 
alpha = params(3); % concentration parameter
%alpha = 1.5; % TODO fix

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

contexts_onehot = zeros(N, D_C);
contexts_onehot(sub2ind([N D_C], 1:N, contexts')) = 1;
contexts = contexts_onehot;


K_C = 0; % # of active context clusters
K_S = 0; % # of active stimulus clusters
% conditionals
P_Zc_given_C = zeros(max_clusters, D_C); % P(Z_c|C) cluster assignments for each context. Note row = C, col = Z
P_Zs_given_S = zeros(max_clusters, D_S); % P(Z_s|S) cluster assignments for each stimulus. Note row = S, col = Z
% marginals
P_Zc = nan(max_clusters); % P(Z_c) = cluster popularities for contexts
P_Zs = nan(max_clusters); % P(Z_s) = cluster popularities for stimuli
% marginals
P_C = nan(max_clusters); % P(C) = context frequencies
P_S = nan(max_clusters); % P(S) = stimulus frequencies
% joint
P_Zc_and_C = nan(max_clusters, D_C); % P(Z_c,C) = P(Z_c|C) P(C) 
P_Zs_and_S = nan(max_clusters, D_S); % P(Z_s,S) = P(Z_s|S) P(S)

Q = Q0 * ones(max_clusters, max_clusters); % Q(C,S) TODO technically V b/c no action

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)

for n = 1:N % for each trial

    r = rewards(n);
    s = find(stimuli(n,:));
    c = find(contexts(n,:));

    if DO_PRINT, fprintf('\n\n\n========== Trial %d: c = %d, s = %d -> r = %d ============\n', n, c, s, r); end

    % Initialize cluster probabilities according to CRP if stimulus and/or context is new
    %
    [P_Zc_given_C, K_C] = CRP_update(P_Zc_given_C, K_C, c, alpha, DO_PRINT);
    [P_Zs_given_S, K_S] = CRP_update(P_Zs_given_S, K_S, s, alpha, DO_PRINT);

    P_C = mean(contexts(1:n,:), 1);
    P_S = mean(stimuli(1:n,:), 1);
    P_Zc_and_C = P_Zc_given_C .* P_C; % ncf matlab won't like this
    P_Zs_and_S = P_Zs_given_S .* P_S;
    P_Zc = sum(P_Zc_and_C, 2);
    P_Zs = sum(P_Zs_and_S, 2);

    if DO_PRINT
        disp('      prior P(Z|C):');
        disp(P_Zc_given_C);
        disp('      prior P(Z|S):');
        disp(P_Zs_given_S);
        disp('      prior Q:');
        disp(Q);
    end
    priors_P_Zc_given_C(:,:,n) = P_Zc_given_C;
    priors_P_Zs_given_S(:,:,n) = P_Zs_given_S;
    prior_Zc_given_c(n,:) = P_Zc_given_C(:,c);
    prior_Zs_given_s(n,:) = P_Zs_given_S(:,s);
    priors_Q(:,:,n) = Q;
    priors_P_Zc_and_C(:,:,n) = P_Zc_and_C;
    priors_P_Zs_and_S(:,:,n) = P_Zs_and_S;
    priors_P_Zc(n,:) = P_Zc;
    priors_P_Zs(n,:) = P_Zs;
    priors_P_C(n,:) = P_C;
    priors_P_S(n,:) = P_S;


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

    % update posterior over clusters for current stimulus and context
    %
    assert(r == 1 || r == 0); % assumes binary rewards
    
    P_Zc_given_C = Reward_update(P_Zc_given_C, K_C, Q(:,z_s), r, c, DO_PRINT);
    P_Zs_given_S = Reward_update(P_Zs_given_S, K_S, Q(z_c,:), r, s, DO_PRINT);

    P_Zc_and_C = P_Zc_given_C .* P_C; % ncf matlab won't like this
    P_Zs_and_S = P_Zs_given_S .* P_S;
    P_Zc = sum(P_Zc_and_C, 2);
    P_Zs = sum(P_Zs_and_S, 2);

    % pick clusters of current stimulus/context for updating (maximum a posteriori)
    %
    [~, z_c] = max(P_Zc_given_C(:,c));
    [~, z_s] = max(P_Zs_given_S(:,s));

    if DO_PRINT, fprintf('   (a posteriori) inferred z_c = %d, z_s = %d, Q = %.4f\n', z_c, z_s, Q(z_c, z_s)); end

    % update Q-value based on reward
    %
    PE = r - Q(z_c, z_s);
    Q(z_c, z_s) = Q(z_c, z_s) + eta * PE;

    if DO_PRINT, fprintf('   PE = %.3f - %.3f = %.3f\n', r, Q(z_c, z_s), PE); end
    if DO_PRINT, fprintf('     new Q(%d, %d) = %.3f\n', z_c, z_s, Q(z_c, z_s)); end

    if DO_PRINT
        disp('      posterior P(Z|C):');
        disp(P_Zc_given_C);
        disp('      posterior P(Z|S):');
        disp(P_Zs_given_S);
        disp('      posterior Q:');
        disp(Q);
    end
    posteriors_Zc_given_C(:,:,n) = P_Zc_given_C;
    posteriors_Zs_given_S(:,:,n) = P_Zs_given_S;
    posterior_Zc_given_c(n,:) = P_Zc_given_C(:,c)';
    posterior_Zs_given_s(n,:) = P_Zs_given_S(:,s)';
    posteriors_Q(:,:,n) = Q;
    PEs(n) = PE;
    posteriors_P_Zc_and_C(:,:,n) = P_Zc_and_C;
    posteriors_P_Zs_and_S(:,:,n) = P_Zs_and_S;
    posteriors_P_Zc(n,:) = P_Zc;
    posteriors_P_Zs(n,:) = P_Zs;
    posteriors_P_C(n,:) = P_C;
    posteriors_P_S(n,:) = P_S;


end


train_results.choices = choices;
train_results.values = values;
train_results.P_C = P_Zc_given_C;
train_results.P_S = P_Zs_given_S;
train_results.K_C = K_C;
train_results.K_S = K_S;
train_results.Q = Q;
train_results.priors_Zc_given_C = priors_P_Zc_given_C;
train_results.priors_Zs_given_S = priors_P_Zs_given_S;
train_results.prior_Zc_given_c = prior_Zc_given_c;
train_results.prior_Zs_given_s = prior_Zs_given_s;
train_results.priors_Q = priors_Q;
train_results.priors_P_Zc_and_C = priors_P_Zc_and_C;
train_results.priors_P_Zs_and_S = priors_P_Zs_and_S;
train_results.priors_P_Zc = priors_P_Zc;
train_results.priors_P_Zs = priors_P_Zs;
train_results.priors_P_C = priors_P_C;
train_results.priors_P_S = priors_P_S;
train_results.posteriors_Zc_given_C = posteriors_Zc_given_C;
train_results.posteriors_Zs_given_S = posteriors_Zs_given_S;
train_results.posterior_Zc_given_c = posterior_Zc_given_c;
train_results.posterior_Zs_given_s = posterior_Zs_given_s;
train_results.posteriors_Q = posteriors_Q;
train_results.PEs = PEs;
train_results.posteriors_P_Zc_and_C = posteriors_P_Zc_and_C;
train_results.posteriors_P_Zs_and_S = posteriors_P_Zs_and_S;
train_results.posteriors_P_Zc = posteriors_P_Zc;
train_results.posteriors_P_Zs = posteriors_P_Zs;
train_results.posteriors_P_C = posteriors_P_C;
train_results.posteriors_P_S = posteriors_P_S;

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

function P = Reward_update(P, K, Q, r, c, DO_PRINT)
    % Update P(Z_c|C) based on reward and Q-values for given c (notice Q-values are a vector here, for a fixed z_s)
    %
    %if DO_PRINT, fprintf('          across context clusters: Q(:, z_s) = %s\n', sprintf('%.3f, ', Q(:,z_s))); end
    for i=1:K % for each context cluster
        q = Q(i);
        lik = binopdf(r, 1, q);
        %if DO_PRINT, fprintf('                 cluster %d: Q(%d, z_s) = %.3f, lik = %.3f, prior P(%d|c) = %.3f\n', i, i, q, lik, i, P_C(c, i)); end
        P(i,c) = P(i,c) * lik;
    end
    P(:,c) = P(:,c) / sum(P(:,c));
end
