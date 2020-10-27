function [lik, out] = ideal_choice(n, particle, stimuli, contexts, rewards, actions, params, which_structures, DO_PRINT)


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

sigma_r = sqrt(0.01);
sigma_w = sqrt(prior_variance); % std for gaussian prior over weights, uncertainty; decreasing the learning rate
tau = sqrt(diffusion_variance);



    % Make prediction for current trial
    % notice conversion to column vectors
    %
    r = rewards(n);
    a = actions(n);
    c = full(ind2vec(contexts(n), K)); % one-hot vector for the context
    x{1} = stimuli(n,:)';
    x{2} = x{1};
    x{3} = [x{1}; c];
    x{4} = c;
    x{5} = c;
    x{6} = 1;

    [V, vals] = model_value_helper(x, particle.w, stimuli(n,:), contexts(n), particle.sample);

    out = predict(V);
    lik = binopdf(a, 1, out);

