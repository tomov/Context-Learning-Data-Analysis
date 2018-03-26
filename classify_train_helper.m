function [classifier, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, subjs, outFilename)

% helper for classify_train.m
%

accuracy = NaN; % TODO make all of them output it
stats = struct;

%
% Fit them
%

switch method
    
    case 'patternnet' % neural network classifier
        
        % patternnet wants column feature vectors. I.e. each data point is a column
        % so we have to rotate it ...
        %
        inputs = inputs'; % ugh MATLAB
        targets = targets';

        % from https://github.com/tomov/food-recognition/blob/master/neural_train.m

        % Create a Pattern Recognition Network
        net = patternnet(10, 'trainscg', 'crossentropy'); % important to use cross-entropy error for multinomial classification

        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 50/100;
        net.divideParam.valRatio = 25/100;
        net.divideParam.testRatio = 25/100;
        net.trainParam.showWindow = false; % don't show GUI on NCF

        % Train the Network
        [net,tr] = train(net,inputs,targets);
        stats.tr = tr;

        % Test the Network
        outputs = net(inputs);

        % Evaluate performance
        errors = gsubtract(targets,outputs);

        %stats.performance = perform(net,targets,outputs);
        %stats.p = 1 - binocdf(size(targets,2) * stats.performance, size(targets,2), 1/size(targets,1));

        %stats.trainPerformance = perform(net, targets(:,tr.trainInd), outputs(:,tr.trainInd));
        %stats.p_train = 1 - binocdf(numel(tr.trainInd) * stats.trainPerformance, numel(tr.trainInd), 1/size(targets,1));

        %stats.valPerformance = perform(net, targets(:,tr.valInd), outputs(:,tr.valInd));
        %stats.p_val = 1 - binocdf(numel(tr.valInd) * stats.valPerformance, numel(tr.valInd), 1/size(targets,1));

        %stats.testPerformance = perform(net, targets(:,tr.testInd), outputs(:,tr.testInd));
        %stats.p_test = 1 - binocdf(numel(tr.testInd) * stats.testPerformance, numel(tr.testInd), 1/size(targets,1));


        % View the Network
        %view(net)

        % View confusion matrix
        %[c,cm,ind,per] = confusion(targets,outputs);

        % patternnet wants column feature vectors, so we rotated those
        % but the rest of the code expects them to be rotated
        % so we undo that...
        % 
        inputs = inputs';
        targets = targets';
        outputs = outputs';

        accuracy = classify_get_accuracy(outputs, targets); % note this include training data too!
        %accuracy = 100 * stats.performance; % WRONG!!!! this is cross-entropy, NOT accuracy!
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
   

    case 'cvpatternnet' % manually cross-validated neural network classifier ; best thing we got so far
        

        [foldid, kfolds] = balanced_folds(runs, subjs, trials, targets);

        inputs = inputs'; % ugh MATLAB
        targets = targets';

        %kfolds = 9;
        %cv = cvpartition(size(inputs,2), 'kfold', kfolds);
        %cv = cvpartition_runs(kfolds, subjs, runs, trials, inputs'); DON'T -- then the folds are unbalanced...

        testOutputs = [];

        for i = 1:kfolds
            % setup network
            net = patternnet(10, 'trainscg', 'crossentropy'); % default params; important to use cross-entropy error for multinomial classification
            net.divideParam.trainRatio = 75/100;
            net.divideParam.valRatio = 25/100;
            net.divideParam.testRatio = 0/100; % we test manually
            net.trainParam.showWindow = false; % don't show GUI on NCF

            % Train the Network
            %train_idx = training(cv,i);
            train_idx = find(foldid ~= i);
            [net,tr] = train(net, inputs(:,train_idx), targets(:,train_idx));
            stats.folds(i).tr = tr;

            % Test the Network
            outputs = net(inputs); % TODO optimize, test only

            %test_idx = test(cv,i); 
            test_idx = find(foldid == i);
            %stats.folds(i).performance = perform(net, targets, outputs);  WARNING -- THIS IS NOT ACCURACY!!! this is the cross-entropy error
            %stats.folds(i).testPerformance = perform(net, targets(:,test_idx), outputs(:,test_idx)); DONT DO IT!
            stats.folds(i).net = net;

            testOutputs(test_idx,:) = outputs(:,test_idx)';
        end

        outputs = testOutputs; % only consider the outputs from the test trials (from all folds)

        accuracy = classify_get_accuracy(outputs, targets');
        k = floor(size(targets,2) * accuracy / 100) - 1; % P(# corr >= ...) = 1 - P(# corr <= ... - 1)
        stats.p = 1 - binocdf(k, size(targets,2), 1/size(targets,1));

        fprintf('Success rate = %.0f%%, p = %f\n', accuracy, stats.p);


        % undo this atrocity
        % 
        inputs = inputs';
        targets = targets';

        classifier = net;



    case 'cvfitcnb' % manually cross-validated naive bayes classifier; for compatibility with searchmight 
        
        [foldid, kfolds] = balanced_folds(runs, subjs, trials, targets);

        testOutputs = [];
        [~, labels] = max(targets, [], 2); % from one-hot vector to indices

        for i = 1:kfolds
            % train
            train_idx = find(foldid ~= i);
            Mdl = fitcnb(inputs(train_idx,:), labels(train_idx));
            
            % test
            [~, outputs] = predict(Mdl, inputs);

            % save outputs
            test_idx = find(foldid == i);
            testOutputs(test_idx,:) = outputs(test_idx,:);
        end

        outputs = testOutputs; % only consider the outputs from the test trials (from all folds)

        accuracy = classify_get_accuracy(outputs, targets);
        k = floor(size(targets,2) * accuracy / 100) - 1; % P(# corr >= ...) = 1 - P(# corr <= ... - 1)
        stats.p = 1 - binocdf(k, size(targets,2), 1/size(targets,1));

        fprintf('Success rate = %.0f%%, p = %f\n', accuracy, stats.p);


        classifier = Mdl;



        
    case 'glmnet' % multinomial GLM classifier
    
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

        % no need to do that manually here -- cvglmnet does it
        % automatically
        %
        %[mses, msesems] = glmnetKFoldCrossValidation(inputs, targets, fitObj, 'multinomial', 'response', 4);
        %[~, lambda_idx] = min(mses); % pick lambda with smallest MSE
        %lambda = fitObj.lambda(lambda_idx);

        outputss = glmnetPredict(fitObj, inputs, fitObj.lambda, 'response');

        accuracies = nan(size(outputss, 3), 1);
        for l = 1:size(outputss, 3) % for all lambdas
            outputs = outputss(:,:,l);
            accuracy = classify_get_accuracy(outputs, targets);
            accuracies(l) = accuracy;
            fprintf('Success rate for %d (lambda = %.4f) = %.2f%%\n', l, fitObj.lambda(l), accuracy);
        end

        classifier = fitObj;
   

        
    case 'cvglmnet' % cross-validated multinomial GLM classifier
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);

        nfolds = 10;

        %c = cvpartition_runs(nfolds, subjs, runs, trials, inputs);
        %foldid = c.Impl.indices;
        %assert(length(foldid) == size(inputs, 1));
        %disp('folds:');
        %disp(foldid);
        
        % x = inputs
        % y = targets
        %
        parallel = false;
        keep = true;
        %CVfit = cvglmnet(inputs, targets, 'multinomial', options, 'deviance', [], foldid, parallel, keep);
        CVfit = cvglmnet(inputs, targets, 'multinomial', options, 'deviance', nfolds, [], parallel, keep);
        %disp(CVfit);

        outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate (lambda = %.4f) is %.2f%%\n', CVfit.lambda_1se, accuracy);

        classifier = CVfit;



    case 'cvcvglmnet' % cross-validated cross-validated multinomial GLM classifier (cvglmnet to pick the lambda's; outer CV for us to test). See cvpatternnet
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);
        nfolds = 5;
        parallel = false;
        keep = true;

        [foldid, kfolds] = balanced_folds(runs, subjs, trials, targets);

        testOutputs = [];
        for i = 1:kfolds
            train_idx = find(foldid ~= i);
            CVfit = cvglmnet(inputs(train_idx,:), targets(train_idx,:), 'multinomial', options, 'deviance', nfolds, [], parallel, keep);

            outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response'); % TODO optimize, test only

            test_idx = find(foldid == i);
            testOutputs(test_idx,:) = outputs(test_idx,:);
        end

        outputs = testOutputs; % only consider the outputs from the test trials (from all folds)

        accuracy = classify_get_accuracy(outputs, targets);
        stats.p = 1 - binocdf(size(targets,2) * accuracy / 100, size(targets,2), 1/size(targets,1));
        fprintf('Success rate = %.0f%%, p = %f\n', accuracy, stats.p);

        classifier = CVfit;



    case 'fitcecoc' % multinomial SVM classifier, cross-validated

        kfolds = 10;
        c = cvpartition_runs(kfolds, subjs, runs, trials, inputs);

        [~, labels] = max(targets, [], 2);
        Mdl = fitcecoc(inputs, labels, 'CVPartition', c, 'FitPosterior', 1, 'Verbose', 2);

        [out_labels, negloss, cost, outputs] = kfoldPredict(Mdl);

        %accuracy = mean(outputs == labels);  lame -- don't use this
        %accuracy = mean(sum(outputs .* targets, 2));
        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate is %.2f%%\n', accuracy);

        classifier = Mdl;


    case 'mnrfit'
        % 
        %
        assert(false);
        
    otherwise
        assert(false, 'should be one of the above');
end


% Save everything except for inputs (too big)
%
if exist('outFilename', 'var') && ~isempty(outFilename)
    fprintf('SAVING to %s\n', outFilename);
    save(outFilename,'-regexp','^(?!(inputs)$).');
end


end



% ensure folds contain whole runs, i.e. no run is split across training & test
% this is to ensure cross-run predictions
% b/c of temporal autocorrelation within runs, it is very easy to classify conditions within runs
%
function c = cvpartition_runs(k, subjs, runs, trials, inputs)
    c_runs = cvpartition(numel(subjs) * numel(runs), 'Kfold', k);
    assert(size(inputs, 1) == numel(subjs) * numel(runs) * numel(trials));
    c = cvpartition(size(inputs, 1), 'Kfold', k);
    i = repmat(c_runs.Impl.indices, 1, numel(trials));
    i = i'; 
    i = i(:);
    c.Impl.indices = i;
    c.Impl.TestSize = accumarray(i, 1)';
    c.Impl.TrainSize = size(i, 1) - c.Impl.TestSize;
end


