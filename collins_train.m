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


K_C = 0; % # of active context clusters
K_S = 0; % # of active stimulus clusters
P_Zc_given_C = zeros(D_C, max_clusters); % (C,Z_c) = P(Z_c|C) cluster assignments for each context. Note row = C, col = Z
P_Zs_given_S = zeros(D_S, max_clusters); % (S,Z_s) = P(Z_s|S) cluster assignments for each stimulus. Note row = S, col = Z
P_C = zeros(max_clusters); % P(Z_c) = cluster popularities for contexts
P_S = zeros(max_clusters); % P(Z_s) = cluster popularities for stimuli
P_C(1) = 1;
P_S(1) = 1;
P_Zc_and_C = zeros(D_C, max_clusters); % P(C,Z_c) = P(Z_c|C) P(C) 
P_Zs_and_S = zeros(D_S, max_clusters); % P(S,Z_s) = P(Z_s|S) P(S)

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
    priors_C(:,:,n) = P_Zc_given_C;
    priors_S(:,:,n) = P_Zs_given_S;
    prior_c(n,:) = P_Zc_given_C(c,:);
    prior_s(n,:) = P_Zs_given_S(s,:);
    priors_Q(:,:,n) = Q;


    % pick clusters of current stimulus/context for action selection (maximum a priori)
    %
    [~, z_c] = max(P_Zc_given_C(c,:));
    [~, z_s] = max(P_Zs_given_S(s,:));

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

    % pick clusters of current stimulus/context for updating (maximum a posteriori)
    %
    [~, z_c] = max(P_Zc_given_C(c,:));
    [~, z_s] = max(P_Zs_given_S(s,:));

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
    posteriors_C(:,:,n) = P_Zc_given_C;
    posteriors_S(:,:,n) = P_Zs_given_S;
    posterior_c(n,:) = P_Zc_given_C(c,:);
    posterior_s(n,:) = P_Zs_given_S(s,:);
    posteriors_Q(:,:,n) = Q;
    PEs(n) = PE;

end


train_results.choices = choices;
train_results.values = values;
train_results.P_C = P_Zc_given_C;
train_results.P_S = P_Zs_given_S;
train_results.K_C = K_C;
train_results.K_S = K_S;
train_results.Q = Q;
train_results.priors_C = priors_C;
train_results.priors_S = priors_S;
train_results.prior_c = prior_c;
train_results.prior_s = prior_s;
train_results.priors_Q = priors_Q;
train_results.posteriors_C = posteriors_C;
train_results.posteriors_S = posteriors_S;
train_results.posterior_c = posterior_c;
train_results.posterior_s = posterior_s;
train_results.posteriors_Q = posteriors_Q;
train_results.PEs = PEs;

end


function [P, K] = CRP_update(P, K, c, alpha, DO_PRINT)
    % if context c is new, update its distribution over clusters P(Z_c|c)
    %
    if sum(P(c,:)) == 0
        % new context
        %
        p = [sum(P(:,1:K)), alpha];
        p = p / sum(p);
        P(c,1:numel(p)) = p;
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
        P(c, i) = P(c, i) * lik;
    end
    P(c, :) = P(c, :) / sum(P(c, :));
end
