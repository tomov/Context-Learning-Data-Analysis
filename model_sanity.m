% Sanity test model_train.m and model_test.m
%

clear;

DO_PLOT = true;

learning_rate = 0.1; % 1 in paper; 0.1 is better
softmax_temp = 2; % 1 in paper; 2 is better

%
% Set up the experiment
%

% irrelevant context group = group 1
% set up stimuli vectors, contexts and rewards for x1c1, x2c1, x2c1, x2c2
%
x{1} = [1 0 0; 0 1 0; 1 0 0; 0 1 0]; % stimuli vectors x1 x2 x1 x2
c{1} = [1; 1; 2; 2]; % context indices c1 c1 c2 c2
r{1} = [1; 0; 1; 0]; % rewards

% modulatory context group = group 2
%
x{2} = x{1};
c{2} = c{1};
r{2} = [1; 0; 0; 1];
% additive context group = group 3
%
x{3} = x{1};
c{3} = c{1};
r{3} = [1; 1; 0; 0];

% repeat trials for each group
%
reps = 5; % = nTrials / 4
for g = 1:3
    x{g} = repmat(x{g}, reps, 1);
    c{g} = repmat(c{g}, reps, 1);
    r{g} = repmat(r{g}, reps, 1);
    assert(size(x{g}, 1) == size(c{g}, 1));
    assert(size(x{g}, 1) == size(r{g}, 1));
end

% test trials
% x1c1, x1c3, x3c1, x3c3
%
test_x = [1 0 0; 1 0 0; 0 0 1; 0 0 1]; % test stimuli: x1 x1 x3 x3
test_c = [1; 3; 1; 3]; % test contexts: c1 c3 c1 c3

Ms = [];

for g=3:3 % for each group
    fprintf('\n\n ---------------- GROUP %d ------------------\n\n', g);

    %[choices, P_n, ww_n, P, ww, values, valuess, likelihoods, new_values, new_valuess, Sigma, lambdas] = 
    train_results = model_train(x{g}, c{g}, r{g}, [learning_rate, softmax_temp], [1 1 1 0], false);

    %[test_choices] =
    test_results = model_test(test_x, test_c, train_results.P_n, train_results.ww_n, [learning_rate, softmax_temp]);
    for n = 1:size(test_x, 1)
        x_n = test_x(n, :)';
        c_n = test_c(n, :);
        out = test_results.choices(n);
        fprintf('TEST predction for x = %d, c = %d is %f\n', find(x_n), c_n, out);
    end
    Ms = [Ms; test_results.choices'];
    
    % Plot posterior probability P(M | h_1:n)
    %
    if DO_PLOT
        figure;

        subplot(3, 2, 1);
        plot(train_results.P, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('P(M | h_{1:n})');
        title('Posterior probability after each trial');
        legend({'M1', 'M2', 'M3', 'M4'});

        % Plot weight matrix ww for M1
        %
        subplot(3, 2, 2);
        plot(train_results.ww_after{1}, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('ww_n');
        title('Weight matrix on each trial for M1');
        legend({'x1', 'x2'});

        % Plot weight matrix ww for M2
        %
        subplot(3, 2, 3);
        plot(train_results.ww_after{2}, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('ww_n');
        title('Weight matrix on each trial for M2');
        legend({'x1c1', 'x2c1', 'c1', 'x1c2', 'x2c2', 'c2'});

        % Plot weight matrix ww for M3
        %
        subplot(3, 2, 4);
        plot(train_results.ww_after{3}, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('ww_n');
        title('Weight matrix on each trial for M3');
        legend({'x1', 'x2', 'c1', 'c2'});    

        % Plot weight matrix ww for M4
        %
        subplot(3, 2, 5);
        plot(train_results.ww_after{4}, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('ww_n');
        title('Weight matrix on each trial for M4');
        legend({'c1', 'c2'});    
    end
end

contexts = {'irrelevant', 'modulatory', 'additive'};

figure;

barweb(Ms, zeros(size(Ms)), 1, contexts, 'Choice probabilities in test phase');
ylabel('Sick probability');
legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'});
