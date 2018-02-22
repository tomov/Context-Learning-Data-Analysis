function train_results = collins_train(stimuli, contexts, rewards, params, which_structures, DO_PRINT)

% Collins & Frank clustering model. Clusters stimuli and contexts independently using DP (CRP).
%

eta = 0.1; % learning rate TODO fit
alpha = 1; % concentration parameter TODO fit
Q0 = 0.5; % prior outcome expectaiton TODO collins = 0.5; Sam = 0; maybe fit?
Q_var = 1; % TODO parameter -- also wtf?
inv_softmax_temp = 2; % beta TODO fit

predict = @(Q_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * Q_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

max_clusters = 4;

D_C = 3;          % # of contexts = K in Sam's model
D_S = size(stimuli, 2); % # of stimuli = D in Sam's model

N = size(stimuli, 1); % # observations


K_C = 0; % # of active context clusters
K_S = 0; % # of active stimulus clusters
P_Z_C = zeros(max_clusters, D_C); % P(Z|C) cluster assignments for each context. Note row = C, col = Z
P_Z_S = zeros(max_clusters, D_S); % P(Z|S) cluster assignments for each stimulus. Note row = S, col = Z

Q = ones(max_clusters, max_clusters); % TODO model with actual actions ? a1, a2 


for i = 1:N

    r = rewards(n);
    s = find(stimuli(n,:));
    c = contexts(n); % one-hot vector for the context

    % Initialize cluster probabilities according to CRP if stimulus and/or context is new
    %
    if sum(P_Z_C(c,:) == 0)
        % new context
        %
        p = [sum(P_Z_C(:,1:K_C)), alpha];
        p = p / sum(p);
        P_Z_C(c,:) = p;
        K_C = K_C + 1; % TODO is this correct?
    end
    if sum(P_Z_S(s,:) == 0)
        % new context
        %
        p = [sum(P_Z_S(:,1:K_S)), alpha];
        p = p / sum(p);
        P_Z_S(s,:) = p;
        K_S = K_S + 1; % TODO is this correct?
    end

    % pick clusters of current stimulus/context (maximum a priori)
    %
    [~, z_c] = max(P_Z_C(c,:));
    [~, z_s] = max(P_Z_S(s,:));

    % pick action
    %
    out = predict(Q(z_c, z_s));
    choices = [choices; out];

    % update Q-value
    %
    PE = r - Q(z_c, z_s);
    Q(z_c, z_s) = Q(z_c, z_s) + eta * PE;

    %
    %
    for i=1:D_C
        lik = normpdf(r, Q(), Q_var);
    end
