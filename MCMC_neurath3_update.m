function particle=ideal_update(n, particle, stimuli, contexts, rewards, params, which_structures, DO_PRINT)


% copy of the loop of model_train.m

assert(numel(params) == 2 || numel(params) == 3 || numel(params) == 4);
prior_variance = params(1);
inv_softmax_temp = params(2);
if numel(params) >= 3
    diffusion_variance = params(3);
else
    diffusion_variance = 0.001;
end
if numel(params) >= 4
    initial_weight = params(4);
else
    initial_weight = 0;
end

predict = @(V_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * V_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

% constants
%
N = size(stimuli, 1); % # of trials
D = size(stimuli, 2); % # of stimuli
K = 3;          % # of contexts
num_structures = 6;

sigma_r = sqrt(0.01);
sigma_w = sqrt(prior_variance); % std for gaussian prior over weights, uncertainty; decreasing the learning rate
tau = sqrt(diffusion_variance);




    % Make prediction for current trial
    % notice conversion to column vectors
    %
    r = rewards(n);
    c = full(ind2vec(contexts(n), K)); % one-hot vector for the context
    x{1} = stimuli(n,:)';
    x{2} = x{1};
    x{3} = [x{1}; c];
    x{4} = c;
    x{5} = c;
    x{6} = 1;
%{ 
    [V, vals] = model_value_helper(x, w, stimuli(n,:), contexts(n), P);

    out = predict(V);
    choices = [choices; out];
    values = [values; V];
    valuess = [valuess; vals];
  %} 
    % log 
    %{
    Sigma_before{1}(:,:,n) = particle.S{1}(1:2,1:2);
    Sigma_before{2}(1:2,1:2,n) = particle.S{2}(1:2,1:2,1);
    Sigma_before{2}(3:4,3:4,n) = particle.S{2}(1:2,1:2,2);
    Sigma_before{3}(:,:,n) = particle.S{3}([1 2 4 5], [1 2 4 5]);
    Sigma_before{4}(:,:,n) = particle.S{4}(1:2,1:2);
    ww_before{1} = [ww_before{1}; particle.w{1}(1:2)'];
    ww_before{2} = [ww_before{2}; reshape(particle.w{2}(1:2,1:2), [1 4])];
    ww_before{3} = [ww_before{3}; particle.w{3}([1:2 4:5])'];
    ww_before{4} = [ww_before{4}; particle.w{4}(1:2)'];
%}
    % Kalman filter update 
    %
    for i = 1:num_structures
        if i == 2
            % for M2, each context has a separate Kalman filter
            % Notice we update all of them but only the active one 'sees' the real stimulus, and we get the prediction from it
            %
            for k = 1:K 
                if k == contexts(n)
                    H = x{i}';
                else
                    H = zeros(1, D); % inactive contexts updated with empty stimulus -> only covariance changes
                end
                F = eye(size(particle.w{i}(:,k), 1));
                B = zeros(size(F));
                u = zeros(size(particle.w{i}(:,k)));
                Q = eye(size(particle.w{i}(:,k), 1)) * tau^2;
                [particle.w{i}(:,k), particle.S{i}(:,:,k), r_mean_temp, r_var_temp] = kalman(r, particle.w{i}(:,k), particle.S{i}(:,:,k), F, B, u, Q, H, sigma_r^2);
                if k == contexts(n) % only consider prediction for current context
                    r_mean{i} = r_mean_temp;
                    r_var{i} = r_var_temp;
                end
            end
        elseif i == 5
            % for M5 = M2', each cue has a separate Kalman filter
            % Notice we update all of them but only the active one 'sees' the real context, and we get the prediction from it
            %
            for d = 1:D
                if stimuli(n,d)
                    H = x{i}';
                else
                    H = zeros(1, K); % inactive cues updated with empty stimulus -> only covariance changes
                end
                F = eye(size(particle.w{i}(:,d), 1));
                B = zeros(size(F));
                u = zeros(size(particle.w{i}(:,d)));
                Q = eye(size(particle.w{i}(:,d), 1)) * tau^2;
                [particle.w{i}(:,d), particle.S{i}(:,:,d), r_mean_temp, r_var_temp] = kalman(r, particle.w{i}(:,d), particle.S{i}(:,:,d), F, B, u, Q, H, sigma_r^2);
                if stimuli(n,d)
                    r_mean{i} = r_mean_temp;
                    r_var{i} = r_var_temp;
                end
            end
        else
            F = eye(size(particle.w{i}, 1));
            B = zeros(size(F));
            u = zeros(size(particle.w{i}));
            Q = eye(size(particle.w{i}, 1)) * tau^2;
            [particle.w{i}, particle.S{i}, r_mean{i}, r_var{i}] = kalman(r, particle.w{i}, particle.S{i}, F, B, u, Q, x{i}', sigma_r^2);
        end 
    end

    % Update posterior over causal structures
    %
    for i = 1:num_structures
        liks(i) = normpdf(r, r_mean{i}, r_var{i});
    end
    particle.P = particle.P .* liks;
    particle.P = particle.P / sum(particle.P);
    if isnan(particle.P(1))
        % edge case -- all likelihoods are 0, or the likelihoods of all
        % non-zero P_n's are 0 => the P_n's become all 0 after this => then
        % we divide by 0 and they all become NaNs
        % in to fix, interpret this scenario as them being all "equally
        % bad" => revert to the uniform
        %
        % TODO SAM ask sam what to do
        %
        particle.P = which_structures / sum(which_structures);
    end
 
    k = poissrnd(1.6); % lambda from Bramley 2017, table 3 exp 1

    % Gibbs sample
    prev_sample = particle.sample;
    for i = 1:k
        ix = find(particle.sample); % current structure
        edge = randsample(2,1); % which edge to sample, 1 = x->r, 2 = c->r
        assert(ismember(ix, [1 2 4 6]));
        if ix == 1
            % M1: x->r
            if edge == 1
                % sample x->r => either M1 or M6
                candidates = [1 0 0 0 0 1];
            else
                % sample c->r => either M1 or M2
                candidates = [1 1 0 0 0 0];
            end 
        elseif ix == 2
            % M2: x->r, c->r
            if edge == 1
                % sample x->r => either M2 or M4
                candidates = [0 1 0 1 0 0];
            else
                % sample c->r => either M1 or M2
                candidates = [1 1 0 0 0 0];
            end 
        elseif ix == 4
            % M4: c->r
            if edge == 1
                % sample x->r => either M2 or M4
                candidates = [0 1 0 1 0 0];
            else
                % sample c->r => either M4 or M6
                candidates = [0 0 0 1 0 1];
            end 
        else
            % M6: none
            if edge == 1
                % sample x->r => either M1 or M6
                candidates = [1 0 0 0 0 1];
            else
                % sample c->r => either M4 or M6
                candidates = [0 0 0 1 0 1];
            end 
        end

        P_prop = particle.P .* candidates; % condition on other edge
        P_prop = P_prop / sum(P_prop); % normalize
        if isnan(P_prop(1))
            % impossible structures => just sample uniformly
            P_prop = candidates / sum(candidates);
        end
        assert(~any(isinf(P_prop)));

        ix = randsample(length(particle.P), 1, true, P_prop);
        particle.sample = zeros(size(particle.P));
        particle.sample(ix) = 1;

        %particle.sample
    end

    if ~all(prev_sample == particle.sample)
        % reset evidence
        particle.P = which_structures / sum(which_structures);
    end



    % log stuff
    %
    %likelihoods = [likelihoods; liks];
    %Posterior = [Posterior; particle.P];
    particle.Posterior = [particle.Posterior; particle.P];
    particle.samples = [particle.samples; particle.sample];