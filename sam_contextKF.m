function results = contextKF(x,r,c,opts)
    
    % Context-learning Kalman filter model.
    %
    % USAGE: results = contextKF(x,r,c,[opts])
    %
    % INPUTS:
    %   x - [N x D] matrix of cues, where each row is a trial and each
    %       column is a cue
    %   r - [N x 1] rewards
    %   c - [N x 1] context indicator, e.g., c(1)=2 means context 2 was
    %       present on trial 1
    %   opts - options structure with the following fields:
    %           .tau2 - transition (diffusion) variance (default: 1)
    %           .sr2 - reward noise variance (default: 0.1)
    %           .sw2 - weight prior variance (default: 1)
    %           .Pm - [1 x 3] prior on models (default: uniform)
    %
    % OUTPUTS:
    %   results - structure containing the following fields:
    %               .Pm - [N x 3] model posterior *before* observing data
    %                     on current trial (i.e., this is the agent's predictive
    %                     knowledge on current trial)
    %               .rhat [N x 1] reward prediction on each trial
    %
    % Sam Gershman, Jan 2016
    
    % parameters
    if nargin < 4 || isempty(opts)
        opts = struct('tau2',1,'sr2',0.1,'sw2',1,'Pm',[1 1 1]/3);
    end
    tau2 = opts.tau2;
    sr2 = opts.sr2;
    sw2 = opts.sw2;
    Pm = opts.Pm;
    
    % initialization
    [N, D] = size(x);
    K = max(c);
    S1 = sw2*eye(D);
    for k=1:K; S2{k} = sw2*eye(D); end
    S3 = sw2*eye(D+K);
    W1 = zeros(D,1);
    W2 = zeros(D,K);
    W3 = zeros(D+K,1);
    
    % run Kalman filter
    for n = 1:N
        
        Pm0 = Pm;
        
        % M1 update
        [W1, S1, lik(1), rhat(1), lambda(1)] = kalman_update(W1,S1,x(n,:)',r(n),tau2,sr2);
        
        % M2 update
        % Run Kalman filter separately for each context. Note that even
        % inactive contexts are updated, but only their covariance changes.
        for k = 1:K
            if c(n) == k; x2 = x(n,:); else x2 = zeros(1,D); end
            [W2(:,k), S2{k}, L, rh, lm] = kalman_update(W2(:,k),S2{k},x2',r(n),tau2,sr2);
            if c(n) == k; lik(2) = L; rhat(2) = rh; lambda(2) = lm; end
        end
        
        % M3 update
        C = zeros(K,1); C(c(n)) = 1; x2 = [x(n,:)'; C]; % augmented stimulus representation
        [W3, S3, lik(3), rhat(3), lambda(3)] = kalman_update(W3,S3,x2,r(n),tau2,sr2);
        
        % update posterior over models
        if all(Pm<1)
            Pm = lik.*Pm;
            if all(Pm==0); Pm = Pm + eps; end
            Pm = Pm./sum(Pm);
        end
        
        % store results
        results.Pm(n,:) = Pm0;
        results.rhat(n,1) = rhat*Pm0';
        results.lambda(n,:) = lambda;
        results.Rhat(n,:) = rhat;
    end
end

%--------------------------%

function [w, S, lik, rhat, lambda] = kalman_update(w,S,x,r,tau2,sr2)
    
    % Kalman filtering equations
    
    D = length(x);
    Q = tau2*eye(D);
    rhat = w'*x;
    S2 = S + Q;
    lambda = x'*S2*x + sr2;
    g = (S2*x)/lambda;
    err = r - rhat;
    w = w + g*err;
    S = S2 - g*x'*S2;
    lik = normpdf(r,rhat,lambda);
end
