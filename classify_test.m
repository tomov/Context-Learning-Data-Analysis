function [inputs, targets, outputs, which_rows] = classify_test(method, classifier, runs, trials, subjs, mask, predict_what, z_score, event)
% Test classifier generated by classify_train.m
%
% INPUT:
% method = 'glmnet', 'patternnet', or 'cvglmnet'
% classifier = object returned by classify_train

rng('shuffle');

fprintf('classify_test\n');
disp(method);

[inputs, targets, which_rows] = classify_get_inputs_and_targets(runs, trials, subjs, mask, predict_what, z_score, event);

[~, maskname, ~] = fileparts(mask);
outFilename = fullfile('classifier', ['classify_test_', method, '_', maskname, '_', predict_what, '_', z_score, '_', random_string(), '.mat']);

switch method
    
    case 'patternnet'
        
        net = classifier;

        % patternnet wants column feature vectors. I.e. each data point is a column
        % so we have to rotate it ...
        %
        inputs = inputs'; % ugh MATLAB
        targets = targets';

        % Test the Network
        outputs = net(inputs);
        errors = gsubtract(targets,outputs);
        performance = perform(net,targets,outputs);

        % View the Network
        %view(net)

        % View confusion matrix
        [c,cm,ind,per] = confusion(targets,outputs);

        % patternnet wants column feature vectors, so we rotated those
        % but the rest of the code expects them to be rotated
        % so we undo that...
        % 
        inputs = inputs';
        targets = targets';
        outputs = outputs';

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate = %.2f%%\n', accuracy);

        % TODO these break for some reason...
        %
        %{
        A = cm';
        B = bsxfun(@rdivide,A,sum(A));
        C = imresize(1 - B, 15, 'method', 'box');
        imshow(C);
        xlabel('Targets');
        ylabel('Outputs');
        %}

    case 'glmnet'
        
        fitObj = classifier;

        outputss = glmnetPredict(fitObj, inputs, fitObj.lambda, 'response');

        for l = 1:size(outputss, 3) % for all lambdas
            outputs = outputss(:,:,l);
            accuracy = classify_get_accuracy(outputs, targets);
            fprintf('Success rate for %d (lambda = %.4f) = %.2f%%\n', l, fitObj.lambda(l), accuracy);
        end

    case 'cvglmnet'
    
        CVfit = classifier;

        outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate (lambda = %.4f) is %.2f%%\n', CVfit.lambda_1se, accuracy);
        
    otherwise
        assert(false, 'should be one of the above');
end

% Save everything except for inputs (too big)
%
%fprintf('SAVING to %s\n', outFilename);
%save(outFilename,'-regexp','^(?!(inputs)$).');
