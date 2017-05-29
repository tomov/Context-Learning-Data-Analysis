function [classifier] = classify_train(method, runs, trials, subjs, mask, predict_what, z_score)
% Train classifier to predict stuff based on neural activity at trial onset
% returns a fitObj that you can pass to glmnetPredict
% or a petternnet net that you can use e.g. like net(inputs)
% for params descriptions, see classify_get_inputs_and_targets
%
% method = 'cvglmnet', 'glmnet', 'patternnet'

rng('shuffle');

fprintf('classify_train\n');
disp(method);

[inputs, targets] = classify_get_inputs_and_targets(runs, trials, subjs, mask, predict_what, z_score);

[~, maskname, ~] = fileparts(mask);
outFilename = fullfile('classifier', ['classify_train_', method, '_', maskname, '_', predict_what, '_', random_string()]);

tic
disp('training classifier...');

%
% Fit them
%

switch method
    
    case 'patternnet'
        
        % patternnet wants column feature vectors. I.e. each data point is a column
        % so we have to rotate it ...
        %
        inputs = inputs'; % ugh MATLAB
        targets = targets';

        % from https://github.com/tomov/food-recognition/blob/master/neural_train.m

        % Create a Pattern Recognition Network
        hiddenLayerSize = 200; % TODO param
        net = patternnet(hiddenLayerSize);

        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;
        net.trainParam.showWindow = false; % don't show GUI on NCF

        % Train the Network
        [net,tr] = train(net,inputs,targets);

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

        classifier = net;

        % TODO these break for some reason...
        %
        %{
        A = cm';
        B = bsxfun(@rdivide,A,sum(A));
        C = imresize(1 - B, 15, 'method', 'box');
        imshow(C);
        xlabel('Targets');
        ylabel('Outputs');


        [b, ib] = sort(diag(B));
        sprintf('%.2f %% bottom 10', mean(b(1:10) * 100))
        sprintf('%.2f %% top 10', mean(b(41:50) * 100))
        sprintf('%.2f %% correct', (1 - c) * 100)
        %}
    
        
    case 'glmnet'
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);

        % x = inputs
        % y = targets
        %
        fitObj = glmnet(inputs, targets, 'multinomial', options);
        glmnetPrint(fitObj);

        %save('classify_glmnet_fitObj_only_1-19_mask_scramble_runs.mat', 'fitObj', 'random_run_labels');

        %[mses, msesems] = glmnetKFoldCrossValidation(inputs, targets, fitObj, 'multinomial', 'response', 4);
        %[~, lambda_idx] = min(mses); % pick lambda with smallest MSE
        %lambda = fitObj.lambda(lambda_idx);

        %save('classify_glmnet_fitObj_only.mat', 'fitObj', 'mses', 'msesems', 'lambda');
        %save('classify_glmnet_w000t.mat', '-v7.3');

        outputss = glmnetPredict(fitObj, inputs, fitObj.lambda, 'response');

        for l = 1:size(outputss, 3) % for all lambdas
            outputs = outputss(:,:,l);
            accuracy = classify_get_accuracy(outputs, targets);
            fprintf('Success rate for %d (lambda = %.4f) = %.2f%%\n', l, fitObj.lambda(l), accuracy);
        end

        classifier = fitObj;
    
        
    case 'cvglmnet'
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);

        % x = inputs
        % y = targets
        %
        CVfit = cvglmnet(inputs, targets, 'multinomial', options, 'deviance');
        disp(CVfit);

        outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate (lambda = %.4f) is %.2f%%\n', CVfit.lambda_1se, accuracy);

        classifier = CVfit;
end

toc

fprintf('SAVING to %s\n', outFilename);
save(outFilename);
