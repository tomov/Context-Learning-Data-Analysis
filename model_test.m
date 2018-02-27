function test_results = model_test(stimuli, contexts, train_results, params)
%function [choices, values, valuess] = model_test(x, k, train_results.P_n, train_results.ww_n, params)

% Kalman filter to predict outcomes based on a learned posterior over
% causal structures and weights for each causal structure. Learning them
% happens in model_train.m
%
% INPUT
% Same as in model_train.m
%
% OUTPUT
% Same as in model_train.m
% 

inv_softmax_temp = params(2);

% constants
%
N = size(stimuli, 1); % # of trials
D = size(stimuli, 2); % # of stimuli
K = 3;          % # of contexts

% predict
predict = @(V_n) 1 ./ (1 + exp(-2 * inv_softmax_temp * V_n + inv_softmax_temp)); % predicts by mapping the expectation to an outcome

choices = []; % history of choices
values = []; % history of predicted outcomes, weighted sum across all models (causal structures)
valuess = []; % history of predicted outcomes, one for each model (causal structure)

w = train_results.ww_n;
P = train_results.P_n;
                    
for n = 1:N
    % Make prediction for current trial
    % notice conversion to column vectors
    %
    c = full(ind2vec(contexts(n), K)); % one-hot vector for the context
    x{1} = stimuli(n,:)';
    x{2} = x{1};
    x{3} = [x{1}; c];
    x{4} = c;
    x{5} = c;
 
    [V, vals] = model_value_helper(x, w, stimuli(n,:), contexts(n), P);

    out = predict(V);
    choices = [choices; out];
    values = [values; V];
    valuess = [valuess; vals];
end

test_results.choices = choices;
test_results.values = values;
test_results.valuess = valuess;

